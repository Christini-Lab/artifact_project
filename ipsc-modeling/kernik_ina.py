import myokit
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


def plot_ap_v():
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
    sim = myokit.Simulation(mod, proto)
    t_max = 8000
    times = np.arange(0, t_max, 1)
    dat = sim.run(t_max, log_times=times)

    axs[0].plot(times, dat['membrane.V'], 'k')
    axs[1].plot(times, dat['ina.i_Na'], 'k')
    axs[1].plot(times, dat['ical.i_CaL'], 'g')
    axs[1].plot(times, dat['ikr.i_Kr'], 'g')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[1].set_xlabel('Time (ms)', fontsize=16)
    axs[0].set_ylabel('Voltage (mV)', fontsize=16)

    plt.show()

plot_ap_v()
