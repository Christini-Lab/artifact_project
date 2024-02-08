from import_rtxi import get_exp_as_df
import myokit
import numpy as np
import matplotlib.pyplot as plt
import xlrd


def get_exp_comp_data(f_name='4_021921_2_alex_control'):
    h5_file = f'./data/h5_files/{f_name}.h5'
    meta = f'./data/meta/{f_name}.xlsx'

    workbook = xlrd.open_workbook(meta)
    sheet = workbook.sheet_by_name('Sheet1')

    trial_names = []
    trial_nums = []

    for rownum in range(27, 37):
        trial = sheet.row_slice(rownum, 0, 4)
        trial_names.append(trial[0].value)
        trial_nums.append(int(trial[1].value))

    min_trial = min(trial_nums)

    for rownum in range(1, 10):
        if sheet.cell(rownum, 6).value == '':
            break
        
        artifact_params = [v.value for v in sheet.row_slice(rownum, 6, 10)]

        if artifact_params[0] > min_trial:
            break
        else:
            cm = artifact_params[1]
            ra = artifact_params[2]
            rm = artifact_params[3]

    exp_dat = {}

    for i, trial_num in enumerate(trial_nums):
        temp_dat = get_exp_as_df(h5_file, trial_num, cm)
        exp_dat[trial_names[i]] = temp_dat['Current (pA/pF)'].values


    return exp_dat, temp_dat['Time (s)'].values, cm, ra, rm


def get_paci_comp_data(cm=None, ra=None, rm=None):
    proto = myokit.load_protocol('./mmt_files/sodium_proto_s_exp.mmt')
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-.080, 10)

    alpha_p_names = ['rspre_0', 'rspre_20', 'rspre_40', 'rspre_60', 'rspre_80']

    paci_dat = {}
    new_params = {'model_parameters.Cm': cm*1E-12, 'voltageclamp.rseries': ra*1E6, 'i_Na.g_Na': 1700}

    for i, alpha_p in enumerate([0, .2, .4, .6, .8]):
        mod = myokit.load_model(f'./mmt_files/paci_artifact.mmt')
        cm = mod['model_parameters']['Cm'].value()
        mod['i_Na']['g_Na_scale'].set_rhs(.3)

        new_params['voltageclamp.alpha_c'] = 0
        new_params['voltageclamp.alpha_p'] = alpha_p

        for k, v in new_params.items():
            if v is None:
                continue
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .0001)

        dat = sim.run(t, log_times=times) #Log only iout 

        i_out = [v/cm for v in dat['voltageclamp.Iout']]

        paci_dat[alpha_p_names[i]] = i_out

    alpha_c_names = ['rscomp_0', 'rscomp_20', 'rscomp_40',
                                'rscomp_60', 'rscomp_80']

    for i, alpha_c in enumerate([0, .2, .4, .6, .8]):
        mod = myokit.load_model(f'./mmt_files/paci_artifact.mmt')

        new_params['voltageclamp.alpha_c'] = alpha_c 
        new_params['voltageclamp.alpha_p'] = .7 

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .0001)

        dat = sim.run(t, log_times=times)

        i_out = [v/cm for v in dat['voltageclamp.Iout']]

        paci_dat[alpha_c_names[i]] = i_out

    return paci_dat, times, dat['voltageclamp.Vc']


def  plot_exp_paci(exp_dat, paci_dat, voltages, which_comp='rscomp_80'):
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    axs[0].plot(voltages)
    axs[1].plot(exp_dat[which_comp])
    axs[1].plot(paci_dat[which_comp])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)
    plt.show()



exp_dat, times_exp, cm, ra, rm = get_exp_comp_data()
paci_dat, times_paci, voltages = get_paci_comp_data(cm=cm,ra=ra)

#plt.plot(times_exp, exp_dat['rscomp_80'])
#plt.plot(times_paci, paci_dat['rscomp_80'])
#plt.show()

plot_exp_paci(exp_dat, paci_dat, voltages)
