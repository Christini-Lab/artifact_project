from deap import base, creator, tools
import myokit
import random
import matplotlib.pyplot as plt
import numpy as np
import pickle
import xlrd
from import_rtxi import get_exp_as_df
from plot_all_models import get_exp_sodium_proto


class MyCustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            creator.create('FitnessMin', base.Fitness, weights=(-1.0,))
            creator.create('Individual', list, fitness=creator.FitnessMin)
            if module == '__main__':
                module = 'vc_opt_ga'
            return super().find_class(module, name)
        except:
            import pdb
            pdb.set_trace()


def get_best_ind(folder):
    with open(f'results/{folder}/pop.pkl', 'rb') as f:
        unpickler = MyCustomUnpickler(f)
        obj = unpickler.load()

    pop = pickle.load(open(f'results/{folder}/pop.pkl', "rb" ))
    last_gen = pop[-1]
    last_gen.sort(key=lambda x: x.fitness.values[0], reverse=True)

    return last_gen[-1]


def simulate_model(ind):
    mod = myokit.load_model('./mmt_files/paci_cardio_lei.mmt')

    for param, val in ind[0].items():
        mod[param.split('.')[0]][param.split('.')[1]].set_rhs(val)

    mod['voltageclamp']['alpha_p'].set_rhs(.75)
    p = mod.get('engine.pace')
    p.set_binding(None)

    #c_m = mod.get('artifact.c_m')
    #c_m.set_rhs(C_M)

    v_cmd = mod.get('voltageclamp.Vc')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-.080, 30)
    sim = myokit.Simulation(mod, holding_proto)
    t_max = holding_proto.characteristic_time()
    sim.run(t_max)
    mod.set_state(sim.state())

    times = np.arange(0, t_max, 0.0001)
    #times = np.arange(0, 2000, 0.1)
    ## CHANGE THIS FROM holding_proto TO SOMETHING ELSE
    #sodium_proto = myokit.load_protocol('./mmt_files/sodium_proto_s_exp.mmt')
    sodium_proto = get_exp_sodium_proto(1)
    t_max = sodium_proto.characteristic_time()
    times = np.arange(0, t_max, 0.0001)
    sim = myokit.Simulation(mod, sodium_proto)

    dat = sim.run(t_max, log_times=times)

    return dat, times


def get_baseline_sim():
    mod, proto, x = myokit.load('./mmt_files/paci.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('Membrane.Vm')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-.080, 30)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    #proto = myokit.load_protocol('./mmt_files/sodium_proto_s_exp.mmt')

    proto = get_exp_sodium_proto(1)
    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .0001)

    dat = sim.run(t, log_times=times)

    return dat, times


def plot_best_mod(folder):
    ind = get_best_ind(folder)
    dat, times = simulate_model(ind)
    dat_base, times_base = get_baseline_sim()
    cm = ind[0]['model_parameters.Cm']
    i_out = [v/cm for v in dat['voltageclamp.Iout']]
    #i_ion_base = [v/28.7E-11 for v in dat_base['Membrane.i_ion']] 

    exp_dat = get_exp_comp_data()
    exp_data = exp_dat[0]['rscomp_80']
    exp_times = exp_dat[1]

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(times, dat['voltageclamp.Vc'])
    axs[1].plot(times_base, dat_base['Membrane.i_ion'], 'k--', label='Baseline')
    axs[1].plot(times, i_out, label='Best fit')
    axs[1].plot(exp_times, exp_data, label='Exp (Rscomp=80, Rspred=70)')

    axs[1].legend()
    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_ylabel('Voltage (V)')
    axs[1].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    print(ind)
    plt.show()

    fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(times, dat['voltageclamp.Vc'])
    axs[1].plot(times_base, dat_base['Membrane.i_ion'], 'k--', label='Baseline')
    axs[1].plot(times, i_out, label='Best fit')
    axs[1].plot(exp_times, exp_data, label='Exp (Rscomp=80, Rspred=70)')
    axs[2].plot(times, dat['i_Na.i_Na'])
    axs[3].plot(times, dat['i_CaL.i_CaL'])
    axs[4].plot(times, dat['i_to.i_to'])

    axs[1].legend()
    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_ylabel('Voltage (V)')
    axs[1].set_ylabel('Current (A/F)')
    axs[2].set_ylabel('INa (A/F)')
    axs[3].set_ylabel('ICaL (A/F)')
    axs[3].set_xlabel('Time (s)')
    print(ind)
    plt.show()


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



plot_best_mod('exp_3')
