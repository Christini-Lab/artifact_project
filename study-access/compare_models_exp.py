import matplotlib.pyplot as plt
import numpy as np
import myokit
from plot_all_models import get_exp_sodium_proto
import xlrd
from import_rtxi import get_exp_as_df
from scipy.signal import find_peaks



#UTILITY FUNCTIONS
def get_exp_sodium_proto(scale=1000):
    proto = myokit.Protocol()
    proto.add_step(-.08*scale, .4*scale)

    for i in range(-70, 25, 10):
        if i == 0:
            proto.add_step(.1/1000*scale, .05*scale)
        else:
            proto.add_step(i/1000*scale, .05*scale)
        proto.add_step(-.08*scale, .4*scale)

    proto.add_step(.04*scale, .05*scale)
    proto.add_step(-.08*scale, .4*scale)

    proto.add_step(.06*scale, .05*scale)
    proto.add_step(-.08*scale, .1*scale)

    return proto


def get_iv_dat(voltages, current):
    iv_dat = {}

    step_idxs = np.where(np.diff(voltages) > .005)[0]

    v_steps = voltages[step_idxs + 10]
    iv_dat['Voltage'] = v_steps

    currs = []
    for idx in step_idxs:
        #start from +3
        temp_currs = current[(idx+3):(idx+103)]
        x = find_peaks(-np.array(temp_currs), distance=5, width=5)
        #if len(x[0]) > 1:

        if len(x[0]) < 1:
            currs.append(np.min(temp_currs))
        else:
            currs.append(temp_currs[x[0][0]])

    iv_dat['Current'] = currs

    return iv_dat


def simulate_model(mod, proto, with_hold=True, sample_freq=0.0001):
    if mod.time_unit().multiplier() == .001:
        scale = 1000
    else:
        scale = 1

    p = mod.get('engine.pace')
    p.set_binding(None)

    v_cmd = mod.get('voltageclamp.Vc')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    if with_hold:
        holding_proto = myokit.Protocol()
        holding_proto.add_step(-.080*scale, 30*scale)
        sim = myokit.Simulation(mod, holding_proto)
        t_max = holding_proto.characteristic_time()
        sim.run(t_max)
        mod.set_state(sim.state())

    t_max = proto.characteristic_time()
    times = np.arange(0, t_max, sample_freq*scale)
    sim = myokit.Simulation(mod, proto)

    dat = sim.run(t_max, log_times=times)

    return dat, times


def get_baseline_iv(model_name, proto, GNa):
    if model_name == 'Paci':
        mod = myokit.load_model('./mmt_files/paci.mmt')
        mem_name = 'Membrane.Vm'
        ion_name = 'Membrane.i_ion'
        scale = 1
        mod['i_Na']['scale_Na'].set_rhs(GNa)
    elif model_name == 'Kernik':
        mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
        mem_name = 'membrane.V'
        ion_name = 'membrane.i_ion'
        scale = 1000 
        mod['ina']['g_scale'].set_rhs(GNa)
    elif model_name == 'Ord':
        mod = myokit.load_model('./mmt_files/ord_na.mmt')
        mem_name = 'Membrane.V'
        ion_name = 'Membrane.i_ion'
        scale = 1000
        mod['INa']['GNa'].set_rhs(GNa)
    else:
        return

    sample_freq = .0001

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get(mem_name)
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-.080*scale, 30*scale)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, sample_freq*scale)

    dat = sim.run(t, log_times=times)

    voltages = np.array(dat[mem_name])
    current = np.array(dat[ion_name]) # in nA

    if model_name == 'Paci':
        voltages *= 1000

    iv_dat = {}

    step_idxs = np.where(np.diff(voltages) > .005)[0]

    v_steps = voltages[step_idxs + 10]
    iv_dat['Voltage'] = v_steps

    currs = []
    for idx in step_idxs:
        temp_currs = current[(idx):(idx+5)]
        x = find_peaks(-np.array(temp_currs), distance=5, width=2)

        if len(x[0]) < 1:
            currs.append(np.min(temp_currs))
        else:
            currs.append(temp_currs[x[0][0]])

    iv_dat['Current'] = currs


    return times, dat, iv_dat


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

    time = temp_dat['Time (s)'].values
    voltages = temp_dat['Voltage (V)'].values

    return exp_dat, time, voltages, cm, ra, rm


def get_current_peak_traces(voltages, current):
    iv_dat = {}

    step_idxs = np.where(np.diff(voltages) > .005)[0]

    v_steps = voltages[step_idxs + 10]

    currs = []
    for idx in step_idxs:
        temp_currs = current[(idx-20):(idx+150)]

        currs.append(temp_currs)

    iv_dat = dict(zip(v_steps, currs))

    return iv_dat


def get_modexp_curr_voltage_dat(f_name, mod_name, comp_set):
    artifact_type = 'lei'
    exp_dat, time, voltages, cm, ra, rm = get_exp_comp_data(f_name)

    currents = exp_dat[comp_set]

    exp_curr_voltage_dat = get_current_peak_traces(voltages, currents)

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    alpha_p, alpha_c = comp_dict[comp_set][0], comp_dict[comp_set][1]

    if mod_name == 'Kernik':
        fi = f'kernik_cardio_{artifact_type}.mmt'
        scale = 1000
        rs = ra * 1E-3
        GNa = .8 
        cm_group = 'geom'
        cm_val = cm 
        gna_group = 'ina.g_scale'
    elif mod_name == 'Ord':
        fi = f'ord_na_{artifact_type}.mmt'
        scale = 1000
        rs = ra * 1E-3 
        GNa = 90 
        cm_group = 'model_parameters'
        cm_val = cm
        gna_group = 'INa.GNa'
    else:
        fi = f'paci_cardio_{artifact_type}.mmt'
        scale = 1
        rs = ra*1E6
        GNa = .6 
        cm_group = 'model_parameters'
        cm_val = cm * 1E-12
        gna_group = 'i_Na.g_Na_scale'

    proto = get_exp_sodium_proto(scale=scale)

    #times, dat = get_baseline_dat(mod_name, proto)

    mod = myokit.load_model(f'./mmt_files/{fi}')
    mod['voltageclamp']['rseries'].set_rhs(rs)
    if 'rseries_est' in mod['voltageclamp']._variables.keys():
        mod['voltageclamp']['rseries_est'].set_rhs(rs)
    mod[cm_group]['Cm'].set_rhs(cm_val)
    mod['voltageclamp']['cm_est'].set_rhs(cm_val)
    mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
    mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)
    mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
    
    dat, times = simulate_model(mod, proto)
    times = times/scale

    if mod_name == 'Paci':
        voltages = np.array(dat['voltageclamp.Vc']) * 1000
        currents = np.array(dat['voltageclamp.Iout']) / cm_val
    else:
        voltages = np.array(dat['voltageclamp.Vc'])
        currents = np.array(dat['voltageclamp.Iout']) / cm_val

    mod_curr_voltage_dat = get_current_peak_traces(voltages, currents)

    return mod_curr_voltage_dat, exp_curr_voltage_dat


#FIGURE FUNCTIONS
def compare_mod_exp(mod_name='Paci', exp_f='4_021921_2_alex_control'):
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    artifact_type = 'lei'
    exp_dat, time, cm, ra, rm = get_exp_comp_data(exp_f)
    axs[1].plot(time, exp_dat['rscomp_80'])
    alpha_p = .7
    alpha_c = .8

    if mod_name == 'Kernik':
        fi = f'kernik_cardio_{artifact_type}.mmt'
        scale = 1000
        rs = ra * 1E-3
        GNa = .8 
        cm_group = 'geom'
        cm_val = cm 
        gna_group = 'ina.g_scale'
    elif mod_name == 'Ord':
        fi = f'ord_na_{artifact_type}.mmt'
        scale = 1000
        rs = ra * 1E-3 
        GNa = 22 
        cm_group = 'model_parameters'
        cm_val = cm
        gna_group = 'INa.GNa'
    else:
        fi = f'paci_cardio_{artifact_type}.mmt'
        scale = 1
        rs = ra*1E6
        GNa = 1
        cm_group = 'model_parameters'
        cm_val = cm * 1E-12
        gna_group = 'i_Na.g_Na_scale'

    proto = get_exp_sodium_proto(scale=scale)

    times, dat = get_baseline_dat(mod_name, proto)

    mod = myokit.load_model(f'./mmt_files/{fi}')
    mod['voltageclamp']['rseries'].set_rhs(rs)
    if 'rseries_est' in mod['voltageclamp']._variables.keys():
        mod['voltageclamp']['rseries_est'].set_rhs(rs)
    mod[cm_group]['Cm'].set_rhs(cm_val)
    mod['voltageclamp']['cm_est'].set_rhs(cm_val)
    mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
    mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)
    mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
    
    dat, times = simulate_model(mod, proto)

    if mod_name == 'Paci':
        voltages = np.array(dat['voltageclamp.Vc']) * 1000
        currents = np.array(dat['voltageclamp.Iout']) / cm_val
    else:
        voltages = np.array(dat['voltageclamp.Vc'])
        currents = np.array(dat['voltageclamp.Iout']) / cm_val
    
    axs[0].plot(times, voltages)
    axs[1].plot(times, currents)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 16
    fig.suptitle(mod_name, fontsize=fs+2)

    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)
    axs[0].set_ylabel('Current (A/F)', fontsize=fs)

    plt.show()


def compare_all_mods_exp(exp_f='4_021921_2_alex_control', comp_set='rscomp_80'):
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    artifact_type = 'lei'
    exp_dat, time, voltages, cm, ra, rm = get_exp_comp_data(exp_f)
    axs[1].plot(time, exp_dat[comp_set], 'k', label=f'Experimental')

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    alpha_p, alpha_c = comp_dict[comp_set][0], comp_dict[comp_set][1]

    models = ['Kernik', 'Ord', 'Paci']

    for mod_name in models:
        if mod_name == 'Kernik':
            fi = f'kernik_cardio_{artifact_type}.mmt'
            scale = 1000
            rs = ra * 1E-3
            GNa = .8 
            cm_group = 'geom'
            cm_val = cm 
            gna_group = 'ina.g_scale'
        elif mod_name == 'Ord':
            fi = f'ord_na_{artifact_type}.mmt'
            scale = 1000
            rs = ra * 1E-3 
            GNa = 90 
            cm_group = 'model_parameters'
            cm_val = cm
            gna_group = 'INa.GNa'
        else:
            fi = f'paci_cardio_{artifact_type}.mmt'
            scale = 1
            rs = ra*1E6
            GNa = .6 
            cm_group = 'model_parameters'
            cm_val = cm * 1E-12
            gna_group = 'i_Na.g_Na_scale'

        proto = get_exp_sodium_proto(scale=scale)

        #times, dat = get_baseline_dat(mod_name, proto)

        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod[cm_group]['Cm'].set_rhs(cm_val)
        mod['voltageclamp']['cm_est'].set_rhs(cm_val)
        mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
        mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)
        mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
        
        dat, times = simulate_model(mod, proto)
        times = times/scale

        if mod_name == 'Paci':
            voltages = np.array(dat['voltageclamp.Vc']) * 1000
            currents = np.array(dat['voltageclamp.Iout']) / cm_val
        else:
            voltages = np.array(dat['voltageclamp.Vc'])
            currents = np.array(dat['voltageclamp.Iout']) / cm_val
        
        axs[1].plot(times, currents, label=mod_name)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].plot(times, voltages, 'k')

    fs = 16
    fig.suptitle(
            f'Models vs Experiment at Pred: {alpha_p}, Comp: {alpha_c}',
            fontsize=fs+2)

    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)
    axs[1].set_ylabel('Current (A/F)', fontsize=fs)

    axs[1].legend()
    plt.show()


def compare_mods_exp_iv(exp_f='4_021921_2_alex_control', comp_set='rscomp_80'):
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    artifact_type = 'lei'
    exp_dat, time, voltages, cm, ra, rm = get_exp_comp_data(exp_f)
    iv_dat = get_iv_dat(voltages, exp_dat[comp_set])
    ax.plot(iv_dat['Voltage']*1000, np.array(iv_dat['Current']),
                                        'k-o', label='Experimental')

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    alpha_p, alpha_c = comp_dict[comp_set][0], comp_dict[comp_set][1]

    models = ['Kernik', 'Ord', 'Paci']

    for mod_name in models:
        if mod_name == 'Kernik':
            fi = f'kernik_cardio_{artifact_type}.mmt'
            scale = 1000
            rs = ra * 1E-3
            GNa = .8 
            cm_group = 'geom'
            cm_val = cm 
            gna_group = 'ina.g_scale'
        elif mod_name == 'Ord':
            fi = f'ord_na_{artifact_type}.mmt'
            scale = 1000
            rs = ra * 1E-3 
            GNa = 90
            cm_group = 'model_parameters'
            cm_val = cm
            gna_group = 'INa.GNa'
        else:
            fi = f'paci_cardio_{artifact_type}.mmt'
            scale = 1
            rs = ra*1E6
            GNa = .6 
            cm_group = 'model_parameters'
            cm_val = cm * 1E-12
            gna_group = 'i_Na.g_Na_scale'

        proto = get_exp_sodium_proto(scale=scale)

        #times, dat = get_baseline_dat(mod_name, proto)

        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod[cm_group]['Cm'].set_rhs(cm_val)
        mod['voltageclamp']['cm_est'].set_rhs(cm_val)
        mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
        mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)
        mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
        
        dat, times = simulate_model(mod, proto)

        if mod_name == 'Paci':
            voltages = np.array(dat['voltageclamp.Vc']) * 1000
            currents = np.array(dat['voltageclamp.Iout']) / cm_val
        else:
            voltages = np.array(dat['voltageclamp.Vc'])
            currents = np.array(dat['voltageclamp.Iout']) / cm_val


        iv_dat = get_iv_dat(voltages, currents)

        ax.plot(iv_dat['Voltage'], np.array(iv_dat['Current']),
                                        '-o', label=mod_name)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    fs = 16
    fig.suptitle(
            f'Models vs Experiment at Pred: {alpha_p}, Comp: {alpha_c}',
            fontsize=fs+2)

    ax.set_ylabel('Current (pA)', fontsize=fs)
    ax.set_xlabel('Voltage (mV)', fontsize=fs)

    ax.legend()
    plt.show()


def compare_all_comp_settings(exp_f='4_021921_2_alex_control'):
    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    fig, axs = plt.subplots(1, 4, figsize=(16, 9), sharey=True)
    artifact_type = 'lei'

    all_vals = list(range(0, 5))
    for k, v in comp_dict.items():
        exp_dat, time, voltages, cm, ra, rm = get_exp_comp_data(exp_f)
        iv_dat = get_iv_dat(voltages, exp_dat[k])
        if 'pre' in k:
            cols = (v[0], v[1], 0)
        else:
            cols = (0, v[0], v[1])
        axs[0].plot(iv_dat['Voltage']*1000, np.array(iv_dat['Current']),
                 '-o', c=cols,
                label=f'Pred: {v[0]}, Comp: {v[1]}')

    axs[0].set_title('Experimental', fontsize=18)

    i = 1
    for mod_name in ['Kernik', 'Ord', 'Paci']:
        if mod_name == 'Kernik':
            fi = f'kernik_cardio_{artifact_type}.mmt'
            scale = 1000
            rs = ra * 1E-3
            GNa = .8 
            cm_group = 'geom'
            cm_val = cm 
            gna_group = 'ina.g_scale'
        elif mod_name == 'Ord':
            fi = f'ord_na_{artifact_type}.mmt'
            scale = 1000
            rs = ra * 1E-3 
            GNa = 110 
            cm_group = 'model_parameters'
            cm_val = cm
            gna_group = 'INa.GNa'
        else:
            fi = f'paci_cardio_{artifact_type}.mmt'
            scale = 1
            rs = ra*1E6
            GNa = .6 
            cm_group = 'model_parameters'
            cm_val = cm * 1E-12
            gna_group = 'i_Na.g_Na_scale'

        proto = get_exp_sodium_proto(scale=scale)

        #times, dat, iv_dat = get_baseline_iv(mod_name, proto, GNa)
        #
        #axs[i].plot(iv_dat['Voltage'], iv_dat['Current'], 'k-o', label='Ideal')
        
        for k, v in comp_dict.items():
            alpha_p, alpha_c = v[0], v[1]
            mod = myokit.load_model(f'./mmt_files/{fi}')
            mod['voltageclamp']['rseries'].set_rhs(rs)
            if 'rseries_est' in mod['voltageclamp']._variables.keys():
                mod['voltageclamp']['rseries_est'].set_rhs(rs)
            mod[cm_group]['Cm'].set_rhs(cm_val)
            mod['voltageclamp']['cm_est'].set_rhs(cm_val)
            mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
            mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)
            mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
            
            dat, times = simulate_model(mod, proto)

            if mod_name == 'Paci':
                voltages = np.array(dat['voltageclamp.Vc']) * 1000
                currents = np.array(dat['voltageclamp.Iout']) / cm_val
            else:
                voltages = np.array(dat['voltageclamp.Vc'])
                currents = np.array(dat['voltageclamp.Iout']) / cm_val

            iv_dat = get_iv_dat(voltages, currents)

            if 'pre' in k:
                cols = (v[0], v[1], 0)
            else:
                cols = (0, v[0], v[1])
            axs[i].plot(iv_dat['Voltage'], np.array(iv_dat['Current']),
                     '-o', c=cols,
                    label=f'Pred: {v[0]}, Comp: {v[1]}')

        axs[i].set_title(mod_name, fontsize=18)

        i+=1
        print(i)


    fs = 16
    for ax in axs:
        ax.set_xlabel('Voltage (mV)', fontsize=fs)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_ylabel('Current (A/F)', fontsize=fs)
    axs[3].legend()

    plt.show()


def plot_multiple_current_traces(f, comp_setting='rscomp_80'):
    # -60, -50, -40, -30, -20, 0, 20, 40, 60
    # For exp and each model:
    # Get the current data from -2 ms to +10 ms after the step to the above voltages
        #

    kernik_dat, exp_dat = get_modexp_curr_voltage_dat(f, 'Kernik', comp_setting)
    ord_dat, exp_dat = get_modexp_curr_voltage_dat(f, 'Ord', comp_setting)
    paci_dat, exp_dat = get_modexp_curr_voltage_dat(f, 'Paci', comp_setting)

    fig, axs = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(12, 9))

    labs = ['Experiment', 'Kernik', 'Ord', 'Paci']
    times = np.linspace(0, 16.9, 170)

    for it, current_dat in enumerate([exp_dat, kernik_dat, ord_dat, paci_dat]):
        i_fig = 0
        for voltage, curr_trace in current_dat.items():
            row = int(i_fig/3)
            col = np.mod(i_fig, 3)

            if voltage in [-70, -.07, -10, -.01, 10, .01]:
                continue
            if it == 0:
                axs[row, col].plot(times, curr_trace, 'k', label=labs[it])
            else:
                axs[row, col].plot(times, curr_trace, label=labs[it])

            if it == 1:
                axs[row, col].set_title(f'Step to: {voltage} mV')

            i_fig += 1


    for row in axs:
        row[0].set_ylabel('Current (A/F)')
        for ax in row:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.yaxis.grid()
            ax.set_ylim(top=15)

    for col in [0, 1, 2]:
        axs[2, col].set_xlabel('Time (ms)')

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    pred_comp = comp_dict[comp_setting]


    fig.suptitle(f'Pred: {int(pred_comp[0]*100)}% & Comp: {int(100*pred_comp[1])}%', fontsize=20)
    axs[2, 2].legend()
    plt.show()


def plot_vest_vcmd():
    ra = 39.3
    cm = 51
    rm = 791.4
    artifact_type = 'lei'
    fi = f'kernik_cardio_{artifact_type}.mmt'
    scale = 1000
    rs = ra * 1E-3
    GNa = .8 
    cm_group = 'geom'
    cm_val = cm 
    gna_group = 'ina.g_scale'

    proto = get_exp_sodium_proto(scale=scale)
    v = [.95, .7]

    alpha_p, alpha_c = v[0], v[1]
    mod = myokit.load_model(f'./mmt_files/{fi}')
    mod['voltageclamp']['rseries'].set_rhs(rs)
    if 'rseries_est' in mod['voltageclamp']._variables.keys():
        mod['voltageclamp']['rseries_est'].set_rhs(rs)
    mod[cm_group]['Cm'].set_rhs(cm_val)
    mod['voltageclamp']['cm_est'].set_rhs(cm_val)
    mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
    mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)
    mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
    
    dat, times = simulate_model(mod, proto)

    fig, axs = plt.subplots(6, 1, sharex=True)

    axs[0].plot(times, dat['voltageclamp.Vc'])
    axs[0].set_ylabel('Vcmd')

    axs[1].plot(times, dat['voltageclamp.Vclamp'])
    axs[1].set_ylabel('Vclamp')

    axs[2].plot(times, dat['voltageclamp.Vest'])
    axs[2].set_ylabel('Vest')

    axs[3].plot(times, dat['membrane.V'])
    axs[3].set_ylabel('Vm')
    
    axs[4].plot(times, dat['voltageclamp.Iout'])
    axs[4].set_ylabel('Iout')

    axs[5].plot(times, dat['voltageclamp.Vc_prime'])
    axs[5].set_ylabel('Vc_prime')

    plt.show()
    

def plot_exp_comp_data(f_name):
    exp_dat, time, voltages, cm, ra, rm = get_exp_comp_data(f_name)

    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axs[0].plot(1000*time, 1000*voltages, 'k')
    col_scale = [.2, .4, .6, .8, 1]

    cols_1 = [[0, v, 0] for v in col_scale]
    cols_2 = [[v, 0, v] for v in col_scale]

    colors = cols_2 + cols_1 

    i = 0
    for k, arr in exp_dat.items():
        axs[1].plot(time*1000, arr, c=colors[i], label=k)
        i+=1 

    axs[0].set_xlim(2649, 2655)
    axs[1].set_ylim(-120, 10)


    plt.show()


    import pdb
    pdb.set_trace()


def main():
    #files with sodium:
    # - '4_021821_1_alex_control' TERRIBLE SEAL - IGNORE
    # - '4_021821_2_alex_control' Pretty good
    # - '4_021821_3_alex_control' BAD SEAL - IGNORE
    # - '4_021821_4_alex_control' MAYBE
    # - '4_021921_1_alex_control' OKAY
    # - '4_021921_2_alex_control' GOOD

    # Function I made to investigate the effects of different artifacts
    #plot_vest_vcmd()


    #compare_mod_exp('Kernik')
    #compare_all_mods_exp(comp_set='rspre_80')
    #compare_mods_exp_iv(comp_set='rspre_80')
    #compare_all_comp_settings('4_021921_2_alex_control')
    #compare_all_comp_settings('4_021821_2_alex_control')
    #compare_all_comp_settings('4_021821_2_alex_control')
    #plot_multiple_current_traces('4_021921_2_alex_control', 'rspre_0')

    #GOOD
    plot_exp_comp_data('4_021821_4_alex_control')
    #
    #plot_exp_comp_data('4_021821_1_alex_control')




if __name__ == "__main__":
    main()

