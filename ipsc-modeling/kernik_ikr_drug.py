import myokit
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

def plot_ap_v():
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    names = ['Baseline', '1nM Dofetilide', '3nM Dofetilide', '6nM Dofetilide']

    for i, block in enumerate([1, .66, .42, .32]):
        mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
        p = mod.get('engine.pace')
        mod['ikr']['g_scale'].set_rhs(block)
        p.set_binding(None)

        v = mod.get('membrane.V')
        v.demote()
        v.set_rhs(0)
        v.set_binding('pace') # Bind to the pacing mechanism

        proto = myokit.Protocol()
        proto.add_step(-80, 1000)
        proto.add_step(40, 1000)
        proto.add_step(-60, 1000)

        t_max = proto.characteristic_time()
        sim = myokit.Simulation(mod, proto)
        times = np.arange(0, t_max, 1)
        dat = sim.run(t_max, log_times=times)

        axs[1].plot(times, dat['ikr.i_Kr'], c=(1-block, 0, 0), label=names[i])

    axs[0].plot(times, dat['membrane.V'], 'k')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=14)



    axs[1].set_xlabel('Time (ms)', fontsize=16)
    axs[0].set_ylabel('Voltage (mV)', fontsize=16)
    axs[1].set_ylabel('Current', fontsize=16)

    plt.legend(fontsize=14)
    plt.show()

def plot_v():
    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))

    names = ['Baseline', '1nM Dofetilide', '3nM Dofetilide', '6nM Dofetilide']

    #for i, block in enumerate([1, .66, .42, .32]):

    for val in [2E-4, 100E-4]:
        mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
        st = mod.state()
        st[2] = val
        mod.set_state(st)

        p = mod.get('engine.pace')
        #mod['ikr']['g_scale'].set_rhs(block)
        #mod['cai']['Cai'].set_rhs(val)
        p.set_binding(None)

        v = mod.get('membrane.V')
        v.demote()
        v.set_rhs(0)
        v.set_binding('pace') # Bind to the pacing mechanism

        proto = myokit.Protocol()
        proto.add_step(-80, 10)
        #proto.add_step(-12, 100)
        #proto.add_step(6, 100)
        #proto.add_step(40, 100)
        proto.add_step(6, 100)
        #proto.add_step(-60, 1000)

        t_max = proto.characteristic_time()
        sim = myokit.Simulation(mod, proto)
        times = np.arange(0, t_max, 1)
        dat = sim.run(t_max, log_times=times)

        #axs[1].plot(times, dat['ikr.i_Kr'], c=(1-block, 0, 0), label=names[i])
        axs[1].plot(times, dat['membrane.i_ion'], label=val)
        axs[0].plot(times, dat['membrane.V'], 'k')
        axs[2].plot(times, dat['cai.Cai'])
        axs[3].plot(times, dat['ical.i_CaL'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=14)



    axs[1].set_xlabel('Time (ms)', fontsize=16)
    axs[0].set_ylabel('Voltage (mV)', fontsize=16)
    axs[1].set_ylabel('Current', fontsize=16)

    axs[1].legend(fontsize=14)
    plt.show()

#plot_ap_v()
plot_v()
