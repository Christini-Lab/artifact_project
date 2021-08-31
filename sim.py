import myokit
import matplotlib.pyplot as plt
import time

import numpy as np

def plot_spont():
    #Plot spontaneous 
    # proto is actually stim, but the length of the stimulation (.001) is so brief that it won't affect the result.
    mod, proto, x = myokit.load('./kernik_mk.mmt')
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(10000)

    fig, axs = plt.subplots(3, 1, True)
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['ikr.IKr'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()


def plot_vc():
    ################
    #Plot VC protocol
    mod = myokit.load_model('./kernik_2019_mc.mmt')
    proto = myokit.load_protocol('./simple_protocol.mmt')
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

    # Add a v1 variable
    v1 = mem.add_variable('v1')
    v1.set_rhs('-150 + 0.1 * engine.time')

    # Add a v2 variable
    v2 = mem.add_variable('v2')
    v2.set_rhs('5694 - 0.4 * engine.time')

    # Add a p variable
    vp = mem.add_variable('vp')
    vp.set_rhs(0)
    vp.set_binding('pace')

    v.set_rhs("""
    piecewise(
        (engine.time >= 300 and engine.time < 700), v1,
        (engine.time >=14410 and engine.time < 14510), v2,
        vp)
    """)

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


def plot_dc():
    #Plot spontaneous 
    # proto is actually stim, but the length of the stimulation (.001) is so brief that it won't affect the result.
    dc_g = .2
    mod, proto, x = myokit.load('./kernik_2019_mc.mmt')
    #proto.pop()
    #proto.schedule(1, 10, 10, 1000, 0)
    mod['parameters']['ik1_ishi_dc_scale'].set_rhs(dc_g)
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(10000)
    
    fig, axs = plt.subplots(3, 1, True)
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['stimulus.i_dc'])


    from cell_models import kernik, protocols

    m = kernik.KernikModel(no_ion_selective_dict={'I_K1_Ishi': dc_g})
    proto = protocols.PacedProtocol('Kernik')
    proto = protocols.SpontaneousProtocol(10000)
    tr = m.generate_response(proto, is_no_ion_selective=True)
    axs[0].plot(tr.t, tr.y)
    axs[1].plot(tr.t, tr.current_response_info.get_current_summed())
    axs[2].plot(tr.t, tr.current_response_info.get_current('I_no_ion'))


    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()


def plot_kernik_mc():
    mod, proto, x = myokit.load('./kernik_2019_mc.mmt')
    sim = myokit.Simulation(mod, proto)
    dat = sim.run(10000)

    fig, axs = plt.subplots(3, 1, True)
    axs[0].plot(dat['engine.time'], dat['membrane.V'])
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    axs[2].plot(dat['engine.time'], dat['ikr.i_Kr'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()


def compare_models():
    mod, proto, x = myokit.load('./kernik_2019_mc.mmt')
    sim = myokit.Simulation(mod, proto)
    dat_1 = sim.run(5000)

    fig, ax = plt.subplots(1, 1, True)

    t = np.loadtxt('./t.csv')
    v = np.loadtxt('./Vm.csv')

    ax.plot(t, v, label='Original-Kernik')
    ax.plot(dat_1['engine.time'], dat_1['membrane.V'], label='Clerx-Kernik')
    #axs[1].plot(dat_1['engine.time'], dat_1['membrane.i_ion'])
    #axs[2].plot(dat_1['engine.time'], dat_1['ikr.i_Kr'])

    #mod, proto, x = myokit.load('./kernik_mk_hs.mmt')
    #sim = myokit.Simulation(mod, proto)
    #dat_2 = sim.run(5000)

    #ax.plot(dat_2['engine.time'], dat_2['membrane.V'], label='Heijman-Kernik')
    #axs[1].plot(dat_2['engine.time'], dat_2['membrane.i_ion'])
    #axs[2].plot(dat_2['engine.time'], dat_2['ikr.IKr'])

    from cell_models import kernik, protocols
    m = kernik.KernikModel()
    p = protocols.SpontaneousProtocol(5000)
    tr = m.generate_response(p, is_no_ion_selective=False)

    ax.plot(tr.t, tr.y, label='Clark-Kernik')
    #axs[1].plot(tr.t, tr.current_response_info.get_current_summed())
    #axs[2].plot(tr.t, tr.current_response_info.get_current('I_Kr'))


    #for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.legend()
    plt.show()


#mod, proto, x = myokit.load('./kernik_2019_mc.mmt')
#sim = myokit.Simulation(mod, proto)
#dat_or = sim.run(50000)
#
#mod, proto, x = myokit.load('./kernik_2019_mc.mmt')
#proto.pop()
#proto.schedule(1, 1000, 10, 1000, 0)
#mod['parameters']['ik1_ishi_dc_scale'].set_rhs(.1)
#sim = myokit.Simulation(mod, proto)
#dat_dc = sim.run(50000)
#
#fig, ax = plt.subplots(1, 1, True)
#
#ax.plot(dat_or['engine.time'], dat_or['membrane.V'], label='No DC')
#ax.plot(dat_dc['engine.time'], dat_dc['membrane.V'], label='With DC')
#
#plt.show()

#plot_dc()
#plot_vc()
#plot_kernik_mc()
compare_models()
