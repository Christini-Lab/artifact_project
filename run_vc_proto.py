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


def prep_model_for_vc(mod):
    # Set up voltage-clamp simulation
    p = mod.get('engine.pace')
    p.set_binding(None)

    # Get membrane potential, demote to an ordinary variable
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
    
    times = np.arange(0, t_max, 0.5)
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

    fig, axs = plt.subplots(7, 1, sharex=True, figsize=(12, 8))
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    #axs[2].plot(dat['engine.time'], dat['ina.i_Na'])

    mod = myokit.load_model('./mmt_files/kernik_artifact.mmt')
    prep_model_for_vc(mod)
    perform_prestep(mod, step_length=30000)
    dat = get_proto_response(mod, './mmt_files/sodium_proto.mmt') 

    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['voltageclamp.Vclamp'])
    axs[3].plot(dat['engine.time'], dat['voltageclamp.Vp'])
    axs[4].plot(dat['engine.time'], dat['voltageclamp.Vest'])
    axs[5].plot(dat['engine.time'], dat['voltageclamp.Iout'])
    axs[6].plot(dat['engine.time'], dat['voltageclamp.Iin'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 14
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)
    axs[2].set_ylabel('Vclamp', fontsize=fs)
    axs[3].set_ylabel('Vp', fontsize=fs)
    axs[4].set_ylabel('Vest', fontsize=fs)
    axs[5].set_ylabel('Iout', fontsize=fs)
    axs[6].set_ylabel('Iin', fontsize=fs)
    axs[6].set_xlabel('Time (ms)', fontsize=fs)

    plt.show()


def main():
    #assess_prestep()
    #run_opt_proto()
    #run_sodium_proto()
    compare_artifact_model()


if __name__ == '__main__':
    main()
