import myokit
import matplotlib.pyplot as plt
import pickle
import time

import numpy as np


def plot_opt_kernik_response():
    mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    all_segments = pickle.load(open('./data/opt_proto_ms.pkl', 'rb'))
    mem = mod.get('membrane')

    piecewise_function = 'piecewise('

    for k, st in all_segments.items():
        v_new = mem.add_variable(k)
        v_new.set_rhs(st[0])

        time_window = f'(engine.time >= {st[1]} and engine.time < {st[2]})'
        piecewise_function += f'({time_window}), {k}, '
        end_time = st[2]

    piecewise_function += 'vp)'
    t_max = end_time

    vp = mem.add_variable('vp')
    vp.set_rhs(0) 

    v = mod.get('membrane.V')

    v.set_binding(None)
    vp.set_binding('pace')

    v.set_rhs(piecewise_function)

    times = np.arange(0, t_max, 0.1)

    sim = myokit.Simulation(mod, proto)

    dat = sim.run(t_max, log_times=times)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    fig.suptitle('Kernik MMT Response', fontsize=18)

    plt.show()


def run_sodium_proto_kernik():
    mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .1)

    dat = sim.run(t, log_times=times)

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
    fig.suptitle('Kernik Baseline', fontsize=18)

    plt.show()


def run_sodium_proto_paci():
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

    proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .0001)

    dat = sim.run(t, log_times=times)

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['Membrane.Vm'])
    axs[1].plot(dat['engine.time'], dat['Membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['i_Na.i_Na'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[2].set_ylabel('INa (pA/pF)', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)
    fig.suptitle('Paci Baseline', fontsize=18)

    plt.show()


def compare_paci_kernik_sodium():
    mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .1)

    dat = sim.run(t, log_times=times)

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], label='Kernik')
    axs[2].plot(dat['engine.time'], dat['ina.i_Na'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)


    #PACI
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

    proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .0001)

    dat = sim.run(t, log_times=times)

    times = [t * 1000 for t in dat['engine.time']]
    v = [v * 1000 for v in dat['Membrane.Vm']]

    axs[1].plot(times, dat['Membrane.i_ion'], label='Paci')
    axs[2].plot(times, dat['i_Na.i_Na'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[2].set_ylabel('INa (pA/pF)', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    axs[1].set_xlim(1980, 2020)

    fig.suptitle('Kernik vs Paci Baseline', fontsize=18)

    axs[1].legend()
    plt.show()


def compare_kernik_artifact_model():
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .1)

    dat = sim.run(t, log_times=times)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], 'k--', label='Ideal')

    for r_series in [.01, .02, .03]:
        mod = myokit.load_model(f'./mmt_files/kernik_artifact.mmt')

        new_params = {'geom.Cm': 20,
                      'voltageclamp.rseries': r_series,
                      'voltageclamp.cprs': 4.188502860812823,
                      'voltageclamp.voffset_eff': 2.9729555590147863,
                      'voltageclamp.gLeak': 0.23719047586907285}

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        holding_proto = myokit.Protocol()
        holding_proto.add_step(-80, 30000)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .1)

        dat = sim.run(t, log_times=times)

        i_out = [v/new_params['geom.Cm'] for v in dat['voltageclamp.Iout']]

        axs[0].plot(dat['engine.time'], dat['membrane.V'])
        axs[1].plot(dat['engine.time'], i_out,
                label=f'Kernik Artifact, Ra={r_series*1000}MOhms')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    axs[1].set_ylim(-150, 40)
    axs[1].set_xlim(1995, 2010)
    axs[1].legend()

    fig.suptitle('Kernik Baseline vs Artifact', fontsize=18)

    plt.show()


def compare_paci_artifact_model(is_normed=True):
    mod = myokit.load_model('./mmt_files/paci.mmt')
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

    proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .0001)

    dat = sim.run(t, log_times=times)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    if not is_normed:
        axs[0].plot(dat['engine.time'], dat['Membrane.Vm'])
        axs[1].plot(dat['engine.time'], dat['Membrane.i_ion'], 'k--', label='Ideal')

    r_vals = [.002, .006, .01, .014, .018, .022]

    for r_series in r_vals:
        mod = myokit.load_model(f'./mmt_files/paci_artifact.mmt')

        new_params = {'voltageclamp.rseries': r_series*1E9}

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        holding_proto = myokit.Protocol()
        holding_proto.add_step(-.080, 30)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .00001)

        dat = sim.run(t, log_times=times)

        i_out = [v/mod['model_parameters']['Cm'].value()
                for v in dat['voltageclamp.Iout']]

        norm_rvals = (r_series - min(r_vals))/max(r_vals)

        if is_normed:
            times = np.array([t for t in dat['engine.time']])
            idx_start = np.argmin(np.abs(times - 1.995))
            idx_end = np.argmin(np.abs(times - 2.010))
            new_times = times[idx_start:idx_end]
            times = times[idx_start:idx_end]
            i_out = np.array(i_out[idx_start:idx_end])

            idx_start = np.argmin(np.abs(times - 2.0002))

            i_out[0:idx_start] = 0
            max_g = np.min(i_out[idx_start:])
            peak_time = times[np.argmin(i_out[idx_start:])]
            
            times -= peak_time

            i_out = i_out / np.abs(max_g)

        #axs[0].plot(times, dat['Membrane.Vm'], c=(0, 0, norm_rvals))
        axs[1].plot(times, i_out,
                label=f'Paci Artifact, Ra={r_series*1000}MOhms',
                c=(norm_rvals, norm_rvals, norm_rvals))

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    if is_normed:
        axs[1].set_xlim(.004, .01)
        axs[1].set_ylim(-1.1, 0)
    else:
        axs[1].set_xlim(1.995, 2.010)
        axs[1].set_ylim(-250, 40)

    axs[1].legend()

    fig.suptitle('Paci Baseline vs Artifact', fontsize=18)

    plt.show()


def compare_paci_conductances_model(is_normed=True):
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))

    g_vals = [.5, 1, 2, 4]

    for g_na in g_vals:
        mod = myokit.load_model(f'./mmt_files/paci_artifact.mmt')

        mod['i_Na']['g_Na'].set_rhs(mod['i_Na']['g_Na'].value() * g_na)

        holding_proto = myokit.Protocol()
        holding_proto.add_step(-.080, 30)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .00005)

        dat = sim.run(t, log_times=times)

        i_out = [v/mod['model_parameters']['Cm'].value()
                for v in dat['voltageclamp.Iout']]

        norm_gvals = (g_na - min(g_vals))/max(g_vals)

        if is_normed:
            times = np.array([t for t in dat['engine.time']])
            idx_start = np.argmin(np.abs(times - 1.995))
            idx_end = np.argmin(np.abs(times - 2.010))
            new_times = times[idx_start:idx_end]
            times = times[idx_start:idx_end]
            i_out = np.array(i_out[idx_start:idx_end])
            voltages = [v for v in dat['Membrane.Vm']]
            voltages = np.array(voltages[idx_start:idx_end])
            i_na = np.array([v for v in dat['i_Na.i_Na']])
            i_na = np.array(i_na[idx_start:idx_end])

            i_na = -i_na / min(i_na)

            idx_start = np.argmin(np.abs(times - 2.0002))

            i_out[0:idx_start] = 0
            max_g = np.min(i_out[idx_start:])
            peak_time = times[np.argmin(i_out[idx_start:])]
            
            times -= peak_time

            i_out = i_out / np.abs(max_g)

        axs[0].plot(times, voltages, c=(norm_gvals, norm_gvals, norm_gvals))
        axs[1].plot(times, i_out,
                label=f'Paci Artifact, G_Na={g_na}',
                c=(norm_gvals, norm_gvals, norm_gvals))
        axs[2].plot(times, i_na, c=(norm_gvals, norm_gvals, norm_gvals))

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[2].set_ylabel('I_Na Current', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    if is_normed:
        axs[1].set_xlim(.004, .01)
        axs[1].set_ylim(-1.1, 0)
    else:
        axs[1].set_xlim(1.995, 2.010)
        axs[1].set_ylim(-250, 40)

    axs[1].legend()

    fig.suptitle('Paci Baseline vs Artifact', fontsize=18)

    plt.show()


def perf_vs_rupture():
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .1)

    dat = sim.run(t, log_times=times)

    iv_ideal = get_iv_dat(dat, is_artifact=False)

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], 'k--', label='Ideal')
    axs[2].plot(dat['engine.time'], dat['ina.i_Na'], 'k--', label='Ideal')

    new_params = {'geom.Cm': 28,
                  'voltageclamp.cprs': 4.188502860812823,
                  'voltageclamp.voffset_eff': -2.8,
                  'voltageclamp.gLeak': 1,
                  'voltageclamp.alpha_c': .7,
                  'voltageclamp.alpha_p': .7}

    cols = ['b', 'g']
    labs = ['Rupture', 'Perforated']
    for i, r_series in enumerate([.005, .03]):
        mod = myokit.load_model(f'./mmt_files/kernik_artifact.mmt')

        new_params['voltageclamp.rseries'] = r_series 
        new_params['voltageclamp.rseries_est'] = r_series * .95

        if labs[i] == 'Rupture':
            new_params['voltageclamp.alpha_c'] = 0
            new_params['voltageclamp.alpha_p'] = 0

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)


        holding_proto = myokit.Protocol()
        holding_proto.add_step(-80, 30000)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .1)

        dat = sim.run(t, log_times=times)

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

    axs[1].set_ylim(-150, 40)
    axs[2].set_ylim(-150, 40)
    axs[1].set_xlim(1995, 2010)
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





def adjust_compensation_kernik(is_exp_dat=True):
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .1)

    dat = sim.run(t, log_times=times)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], 'k--', label='Ideal')

    new_params = {}

    for alpha_p in [0, .2, .4, .6, .8]:
        mod = myokit.load_model(f'./mmt_files/kernik_artifact.mmt')

        new_params['voltageclamp.alpha_c'] = 0
        new_params['voltageclamp.alpha_p'] = alpha_p

        for k, v in new_params.items():
            par_chil = k.split('.')
            mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        t = holding_proto.characteristic_time()
        sim = myokit.Simulation(mod, holding_proto)
        dat = sim.run(t)

        mod.set_state(sim.state())

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()
        times = np.arange(0, t, .1)

        dat = sim.run(t, log_times=times)


        i_out = [v/mod['geom']['Cm'].value() for v in dat['voltageclamp.Iout']]

        num = alpha_p + .2
        axs[1].plot(dat['engine.time'], i_out, label=f'pred={alpha_p}',
                c=(alpha_p, 0, alpha_p))


    for alpha_c in [0, .2, .4, .6, .8]:
        mod = myokit.load_model(f'./mmt_files/kernik_artifact.mmt')

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
        times = np.arange(0, t, .1)

        dat = sim.run(t, log_times=times)

        i_out = [v/mod['geom']['Cm'].value() for v in dat['voltageclamp.Iout']]

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


def adjust_compensation_paci(is_exp_dat=True):
    mod = myokit.load_model('./mmt_files/paci.mmt')
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

    proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .0001)

    dat = sim.run(t, log_times=times)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['Membrane.Vm'])
    axs[1].plot(dat['engine.time'], dat['Membrane.i_ion'], 'k--', label='Ideal')

    new_params = {}

    cm = mod['model_parameters']['Cm'].value()
    for alpha_p in [0, .2, .4, .6, .8]:
        mod = myokit.load_model(f'./mmt_files/paci_artifact.mmt')

        new_params['voltageclamp.alpha_c'] = 0
        new_params['voltageclamp.alpha_p'] = alpha_p

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

        num = alpha_p + .2
        axs[1].plot(dat['engine.time'], i_out, label=f'pred={alpha_p}',
                c=(alpha_p, 0, alpha_p))


    for alpha_c in [0, .2, .4, .6, .8]:
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

        axs[1].plot(dat['engine.time'], i_out, label=f'comp={alpha_c}',
                c=(0, alpha_c, 0))
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


def compare_paci_kernik():
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))

    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    proto = myokit.load_protocol('./mmt_files/sodium_proto.mmt')

    vm_names = ['membrane.V', 'Membrane.Vm']
    t_step = [.05, .0001]
    i_ion_name = ['membrane.i_ion', 'Membrane.i_ion']

    for i, mod_proto in enumerate([['kernik_2019_mc.mmt', 'sodium_proto.mmt'],
                                    ['paci.mmt', 'sodium_proto_s.mmt']]):
        mod = myokit.load_model(f'./mmt_files/{mod_proto[0]}')
        proto = myokit.load_protocol(f'./mmt_files/{mod_proto[1]}')

        p = mod.get('engine.pace')
        p.set_binding(None)

        # Get membrane potential, demote to an ordinary variable
        #if is_v_demoted:
        v = mod.get(vm_names[i])
        v.demote()
        v.set_rhs(0)
        v.set_binding('pace') # Bind to the pacing mechanism

        sim = myokit.Simulation(mod, proto)
        t_max = proto.characteristic_time()
        times = np.arange(0, t_max, t_step[i])

        dat = sim.run(t_max, log_times=times)

        if i > 0:
            times = times * 1000
        axs[1+i].plot(times,
                dat[i_ion_name[i]], 'k--', label='Ideal')

    v = [v*1000 for v in dat['Membrane.Vm']]
    axs[0].plot(times, v)

    for alpha_p in [0, .2, .4, .6, .8]:
        for i, mod_proto in enumerate([['kernik_artifact.mmt', 'sodium_proto.mmt'],
                            ['paci_artifact.mmt', 'sodium_proto_s.mmt']]):
            mod = myokit.load_model(f'./mmt_files/{mod_proto[0]}')
            proto = myokit.load_protocol(f'./mmt_files/{mod_proto[1]}')

            mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
            mod['voltageclamp']['alpha_c'].set_rhs(0)

            sim = myokit.Simulation(mod, proto)
            t_max = proto.characteristic_time()
            times = np.arange(0, t_max, t_step[i])

            dat = sim.run(t_max, log_times=times)

            if i > 0:
                times = times * 1000
                cm = 2.8E-11
                i_out = [i_out/cm for i_out in dat['voltageclamp.Iout']]
            else:
                cm = 28
                i_out = [i_out/cm for i_out in dat['voltageclamp.Iout']]

            axs[1+i].plot(times,
                    i_out, label=f'Pred={alpha_p*100}',
                    c=(alpha_p, 0, alpha_p))


    for alpha_c in [0, .2, .4, .6, .8]:
        for i, mod_proto in enumerate([['kernik_artifact.mmt', 'sodium_proto.mmt'],
                            ['paci_artifact.mmt', 'sodium_proto_s.mmt']]):
            mod = myokit.load_model(f'./mmt_files/{mod_proto[0]}')
            proto = myokit.load_protocol(f'./mmt_files/{mod_proto[1]}')

            mod['voltageclamp']['alpha_p'].set_rhs(.7)
            mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)

            sim = myokit.Simulation(mod, proto)
            t_max = proto.characteristic_time()
            times = np.arange(0, t_max, t_step[i])

            dat = sim.run(t_max, log_times=times)

            if i > 0:
                times = times * 1000
                cm = 2.8E-11
                i_out = [i_out/cm for i_out in dat['voltageclamp.Iout']]
            else:
                cm = 28
                i_out = [i_out/cm for i_out in dat['voltageclamp.Iout']]

            axs[1+i].plot(times,
                    i_out, label=f'Pred={alpha_c*100}',
                    c=(0, alpha_c, 0))

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel("Command Voltage", fontsize=fs)
    axs[1].set_ylabel("Kernik Iout", fontsize=fs)
    axs[2].set_ylabel("Paci Iout", fontsize=fs)
    axs[2].set_xlabel("Time (ms)", fontsize=fs)

    axs[1].set_ylim(-200, 10)
    axs[2].set_ylim(-200, 10)

    axs[2].legend()
    plt.show()


def main():
    #plot_opt_kernik_response()
    run_sodium_proto_kernik()
    #run_sodium_proto_paci()
    #compare_paci_kernik_sodium()
    #compare_kernik_artifact_model()
    #compare_paci_artifact_model()
    #compare_paci_conductances_model()
    #adjust_compensation_kernik()
    #adjust_compensation_paci()



if __name__ == '__main__':
    main()
