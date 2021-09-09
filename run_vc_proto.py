import myokit
import matplotlib.pyplot as plt
import pickle
import time

import numpy as np


def get_opt_response(mod, proto):
    mem = mod.get('membrane')
    all_segments = pickle.load(open('./data/proto.pkl', 'rb'))
    piecewise_function = 'piecewise('

    for k, st in all_segments.items():
        v_new = mem.add_variable(k)
        v_new.set_rhs(st[0])

        time_window = f'(engine.time >= {st[1]} and engine.time < {st[2]})'
        piecewise_function += f'({time_window}), {k}, '
        end_time = st[2]

    piecewise_function += 'vp)'
    t_max = end_time

    # Add a p variable
    vp = mem.add_variable('vp')
    vp.set_rhs(0) 

    v = mod.get('membrane.V')
    v.set_binding(None)
    vp.set_binding('pace')

    v.set_rhs(piecewise_function)

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)

    dat = sim.run(t_max, log_times=times)

    return dat


def perform_prestep(mod, step_length=30000):
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, step_length)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())


def prep_model_for_vc(mod, is_v_demoted=True):
    # Set up voltage-clamp simulation
    p = mod.get('engine.pace')
    p.set_binding(None)

    # Get membrane potential, demote to an ordinary variable
    #if is_v_demoted:
    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism


def assess_prestep():
    # Load model and pre-step proto
    mod1, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
    mod2, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')

    prep_model_for_vc(mod1)
    prep_model_for_vc(mod2)

    #perform_prestep(mod1)

    dat1 = get_opt_response(mod1, proto)
    dat2 = get_opt_response(mod2, proto)

    fig, axs = plt.subplots(3, 1, True)
    axs[0].plot(dat1['engine.time'], dat1['membrane.V'])
    axs[1].plot(dat1['engine.time'], dat1['membrane.i_ion'])
    axs[2].plot(dat1['engine.time'], dat1['ikr.i_Kr'])

    axs[0].plot(dat2['engine.time'], dat2['membrane.V'])
    axs[1].plot(dat2['engine.time'], dat2['membrane.i_ion'])
    axs[2].plot(dat2['engine.time'], dat2['ikr.i_Kr'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()


def run_opt_proto():
    mod1, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
    prep_model_for_vc(mod1)
    perform_prestep(mod1, step_length=30000)

    dat1 = get_opt_response(mod1, proto)
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat1['engine.time'], dat1['membrane.V'])
    axs[1].plot(dat1['engine.time'], dat1['membrane.i_ion'])
    axs[2].plot(dat1['engine.time'], dat1['ikr.i_Kr'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[2].set_ylabel('IKr (pA/pF)', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    plt.show()


def get_proto_response(mod, proto_location):
    proto = myokit.load_protocol(proto_location)
    t_max = proto.characteristic_time()
    sim = myokit.Simulation(mod, proto)
    
    times = np.arange(0, t_max, 0.05)
    dat = sim.run(t_max, log_times=times)

    return dat


def run_sodium_proto():
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    prep_model_for_vc(mod)
    perform_prestep(mod, step_length=30000)
    dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['ina.i_Na'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[2].set_ylabel('INa (pA/pF)', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    plt.show()


def compare_artifact_model():
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    prep_model_for_vc(mod)
    perform_prestep(mod, step_length=30000)
    dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], label='Ideal')

    for f in ['kernik_artifact.mmt', 'kernik_artifact2.mmt']:
        mod = myokit.load_model(f'./mmt_files/{f}')


        new_params = {'geom.Cm': 20,
                      'voltageclamp.rseries': 0.012054223352242087,
                      'voltageclamp.cprs': 4.188502860812823,
                      'voltageclamp.voffset_eff': 2.9729555590147863,
                      'voltageclamp.gLeak': 0.23719047586907285}

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        #prep_model_for_vc(mod, False)
        perform_prestep(mod, step_length=30000)
        dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

        i_in = [v/new_params['geom.Cm'] for v in dat['voltageclamp.Iin']]
        i_out = [v/new_params['geom.Cm'] for v in dat['voltageclamp.Iout']]

        axs[0].plot(dat['engine.time'], dat['membrane.V'])
        axs[1].plot(dat['engine.time'], i_out, label=f)
        #axs[1].plot(dat['engine.time'], i_in)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    axs[1].set_ylim(-200, 40)
    axs[1].legend()

    plt.show()


def adjust_compensation(is_exp_dat=True):
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    prep_model_for_vc(mod)
    perform_prestep(mod, step_length=30000)
    dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], 'k--', label='Ideal')

    new_params = {'geom.Cm': 28,
                  'voltageclamp.rseries': 0.020,
                  'voltageclamp.rseries_est': 0.019,
                  'voltageclamp.cprs': 4.188502860812823,
                  'voltageclamp.voffset_eff': 2.9729555590147863,
                  'voltageclamp.gLeak': 0.23719047586907285}

    for alpha_p in [0, .2, .4, .6, .8]:
        mod = myokit.load_model(f'./mmt_files/kernik_artifact2.mmt')

        new_params['voltageclamp.alpha_c'] = 0
        new_params['voltageclamp.alpha_p'] = alpha_p

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        perform_prestep(mod, step_length=30000)
        dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

        i_out = [v/new_params['geom.Cm'] for v in dat['voltageclamp.Iout']]

        num = alpha_p + .2
        axs[1].plot(dat['engine.time'], i_out, label=f'pred={alpha_p}',
                c=(alpha_p, 0, alpha_p))


    for alpha_c in [0, .2, .4, .6, .8]:
        mod = myokit.load_model(f'./mmt_files/kernik_artifact2.mmt')

        new_params['voltageclamp.alpha_c'] = alpha_c 
        new_params['voltageclamp.alpha_p'] = .7 

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        #prep_model_for_vc(mod, False)
        perform_prestep(mod, step_length=30000)
        dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

        i_out = [v/new_params['geom.Cm'] for v in dat['voltageclamp.Iout']]

        axs[1].plot(dat['engine.time'], i_out, label=f'comp={alpha_c}',
                c=(0, alpha_c, 0))
        #axs[1].plot(dat['engine.time'], i_in)

    if is_exp_dat:
        '4_021821_4_alex_control.xlsx'
        #plot exp data

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    axs[1].set_ylim(-200, 40)
    axs[1].legend()

    plt.show()


def perf_vs_rupture():
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    prep_model_for_vc(mod)
    perform_prestep(mod, step_length=30000)
    dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 
    iv_ideal = get_iv_dat(dat, is_artifact=False)

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], 'k--', label='Ideal')
    axs[2].plot(dat['engine.time'], dat['ina.i_Na'], 'k--', label='Ideal')

    new_params = {'geom.Cm': 28,
                  'voltageclamp.cprs': 4.188502860812823,
                  'voltageclamp.voffset_eff': -2.8,
                  'voltageclamp.gLeak': 1,
                  'voltageclamp.alpha_c': 0,
                  'voltageclamp.alpha_p': 0}

    cols = ['b', 'g']
    labs = ['Rupture', 'Perforated']
    for i, r_series in enumerate([.003, .03]):
        mod = myokit.load_model(f'./mmt_files/kernik_artifact2.mmt')

        new_params['voltageclamp.rseries'] = r_series 
        new_params['voltageclamp.rseries_est'] = r_series * .95

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        perform_prestep(mod, step_length=30000)
        dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

        i_out = [v/new_params['geom.Cm'] for v in dat['voltageclamp.Iout']]

        axs[0].plot(dat['engine.time'], dat['membrane.V'], c=cols[i])
        axs[1].plot(dat['engine.time'], i_out, c=cols[i],
                label=f'Rseries={r_series*1000}M – {labs[i]}')
        axs[2].plot(dat['engine.time'], dat['ina.i_Na'], c=cols[i])
        if i == 0:
            iv_rupt = get_iv_dat(dat, is_artifact=True)
        else:
            iv_perf = get_iv_dat(dat, is_artifact=True)

    #if is_exp_dat:
    #    '4_021821_4_alex_control.xlsx'
    #    #plot exp data

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[2].set_ylabel('INa', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    axs[1].set_ylim(-200, 40)
    axs[1].legend()

    plt.show()


    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    labs = ['Ideal', 'Rupture', 'Perforated']
    cols = ['k', 'b', 'g']

    for i, iv_dat in enumerate([iv_ideal, iv_rupt, iv_perf]):
        ax.plot(iv_dat[0], iv_dat[1], '-o', c=cols[i], label=labs[i])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    fs = 14
    ax.set_xlabel('Voltage', fontsize=fs)
    ax.set_ylabel('Current', fontsize=fs)
        
    plt.legend()
    plt.show()


def get_iv_dat(dat, is_artifact):
    cm = 28

    times = dat['engine.time']

    if is_artifact:
        i_out = [v / cm for v in dat['voltageclamp.Iout']]
        voltage = dat['voltageclamp.Vc']
    else:
        i_out = dat['membrane.i_ion']
        voltage = dat['membrane.V']

    freq = 1 / (times[1] - times[0])
    voltages = []
    peak_currents = []

    for time in np.linspace(400, 5600, 14):
        st_idx = freq * time
        v = voltage[int(st_idx + 3)]
        voltages.append(v)

        if not is_artifact:
            max_i = np.min(i_out[int(st_idx):int(st_idx+freq*5)])
        else:
            max_i = np.min(i_out[int(st_idx+.2*freq):int(st_idx+freq*5)])
        peak_currents.append(max_i)

    return [voltages, peak_currents] 


def main():
    #assess_prestep()
    #run_opt_proto()
    #run_sodium_proto()
    #compare_artifact_model()
    adjust_compensation()
    #perf_vs_rupture()


if __name__ == '__main__':
    main()
