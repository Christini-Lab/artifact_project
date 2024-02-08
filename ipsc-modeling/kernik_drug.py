import myokit
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks



def plot_v_ca():
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12, 8))
    mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
    sim = myokit.Simulation(mod, proto)
    t_max = 2000
    times = np.arange(0, t_max, 1)
    dat = sim.run(t_max, log_times=times)

    v = np.array(dat['membrane.V'])
    v = (v - v.min()) / (v.max() - v.min())

    ca = np.array(dat['cai.Cai'])
    ca = (ca - ca.min()) / (ca.max() - .3*ca.min())

    ax.plot(times, v, 'k')
    ax.plot(times, ca, 'y--')
    ax.set_xlim(200, 1350)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis=u'y', which=u'both',length=0)
    ax.tick_params(labelleft=False)    


    plt.show()


def plot_ap_v():
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12, 8))
    mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
    sim = myokit.Simulation(mod, proto)
    t_max = 8000
    times = np.arange(0, t_max, 1)
    dat = sim.run(t_max, log_times=times)

    ap_t, ap_v = return_ap(times, dat)
    ax.plot(ap_t, ap_v, 'k')
    ax.set_xlim(-100, 600)
    ax.set_ylim(-80, 40)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Time (ms)', fontsize=16)
    ax.set_ylabel('Voltage (mV)', fontsize=16)

    ax.tick_params(axis=u'both', labelsize=14)

    plt.savefig('figs/baseline_ap.png')
    plt.show()


def plot_ap_ikr_block():
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12, 8))
    st = ['k', 'r--']

    for i, block in enumerate([1, .8]):
        mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
        mod['ikr']['g_scale'].set_rhs(block)
        sim = myokit.Simulation(mod, proto)
        t_max = 8000
        times = np.arange(0, t_max, 1)
        dat = sim.run(t_max, log_times=times)

        ap_t, ap_v = return_ap(times, dat)
        ax.plot(ap_t, ap_v, st[i])

    ax.set_xlim(-100, 600)
    ax.set_ylim(-80, 40)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Time (ms)', fontsize=16)
    ax.set_ylabel('Voltage (mV)', fontsize=16)

    ax.tick_params(axis=u'both', labelsize=14)

    plt.savefig('./figs/ikr_ap_block.png')
    plt.show()


def plot_ap_ical_block():
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12, 8))
    st = ['k', 'r--', 'b--']

    for i, block in enumerate([1, .8, .8]):
        mod, proto, x = myokit.load('./mmt_files/kernik_2019_mc.mmt')
        mod['ikr']['g_scale'].set_rhs(block)
        if i == 2:
            mod['ical']['g_scale'].set_rhs(.2)

        sim = myokit.Simulation(mod, proto)
        t_max = 8000
        times = np.arange(0, t_max, 1)
        dat = sim.run(t_max, log_times=times)

        ap_t, ap_v = return_ap(times, dat)

        ax.plot(ap_t, ap_v, st[i])

    ax.set_xlim(-100, 600)
    ax.set_ylim(-80, 40)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel('Time (ms)', fontsize=16)
    ax.set_ylabel('Voltage (mV)', fontsize=16)

    ax.tick_params(axis=u'both', labelsize=14)

    plt.savefig('./figs/ikr_ical_ap_block.png')
    plt.show()


def return_ap(times, dat):
    v = np.array(dat['membrane.V'])
    pks = find_peaks(np.diff(v), 4)

    start = times[pks[0][-3]] - 200
    cycle = times[pks[0][-2]] - times[pks[0][-3]]
    end = start + cycle

    start_idx = np.argmin(np.abs(times - start))
    end_idx = np.argmin(np.abs(times - end))

    t = times[start_idx:end_idx] - pks[0][-3]
    v = v[start_idx:end_idx]

    return t, v

#ical
plot_ap_v()
plot_ap_ikr_block()
plot_ap_ical_block()
