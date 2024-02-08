import myokit
import matplotlib.pyplot as plt
import time

import numpy as np


def plot_vc():
    ################
    #Plot VC protocol
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
    proto = myokit.load_protocol('./mmt_files/simple_protocol.mmt')
    proto = myokit.load_protocol('./mmt_files/holding_proto.mmt')
    t_max = proto.characteristic_time()
    mem = mod.get('membrane')

    #Set up voltage-clamp simulation
    # Get pacing variable, remove binding
    p = mod.get('engine.pace')
    p.set_binding(None)
    # Get membrane potential, demote to an ordinary variable
    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    #v.set_binding('pace') # Bind to the pacing mechanism

    piecewise_function = 'piecewise('
    
    for i in range(0, 4):
        v_new = mem.add_variable(f'v{i}')
        if (i % 2 == 0):
            v_new.set_rhs(f'{np.random.normal()*10} + engine.time * {.1*np.random.normal()}')
        else:
            v_new.set_rhs(f'{np.random.normal()*10}')

        
        piecewise_function += f'(engine.time >= {i*1000} and engine.time < {(i + 1) * 1000}), v{i}, '

    # Add a v1 variable
    #v1 = mem.add_variable('v1')
    #v1.set_rhs('-150 + 0.1 * engine.time')

    ## Add a v2 variable
    #v2 = mem.add_variable('v2')
    #v2.set_rhs('5694 - 0.4 * engine.time')


    piecewise_function += 'vp)'

    # Add a p variable
    vp = mem.add_variable('vp')
    vp.set_rhs(0)
    vp.set_binding('pace')

    v.set_rhs(piecewise_function)
    #v.set_rhs("""
    #piecewise(
    #    (engine.time >= 300 and engine.time < 700), v1,
    #    (engine.time >=14410 and engine.time < 14510), v2,
    #    vp)
    #""")

    times = np.arange(0, t_max, 0.1)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(t_max, log_times=times)

    fig, axs = plt.subplots(3, 1, True)
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['ikr.i_Kr'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()

plot_vc()
