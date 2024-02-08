import myokit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def get_kernik_response(proto, g_leak):
    mod = myokit.load_model('./mmt_files/kernik_leak.mmt')

    mod['membrane']['gLeak'].set_rhs(g_leak)

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()

    times = np.arange(0, t, 1)

    dat = sim.run(t, log_times=times)

    return dat

    #fig, axs = plt.subplots(2, 1, sharex=True)

    #axs[0].plot(dat['engine.time'], dat['membrane.V'])
    #axs[1].plot(dat['engine.time'], dat['membrane.i_ion'])
    #plt.show()

proto = myokit.Protocol()

for voltage in range(-120, 60, 10):
    if voltage == 0:
        proto.add_step(voltage+1, 2000)
    else:
        proto.add_step(voltage, 2000)

#plot_proto_response(proto)


dat_50 = get_kernik_response(proto, .5)
dat_100 = get_kernik_response(proto, 1)
dat_200 = get_kernik_response(proto, 2)

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))


labs = [r'$g_{leak}$=.5G',r'$g_{leak}$=1G',r'$g_{leak}$=2G']
all_iv_currs = []
for i, dat in enumerate([dat_200, dat_100, dat_50]):
    axs[1].plot(dat['engine.time'], dat['membrane.i_ion'], label=labs[i])
    iv_curr = []
    for i, v in enumerate(list(range(-120, 60, 10))):
        curr_t = i * 2000 +1900
        iv_curr.append(dat['membrane.i_ion'][int(curr_t)])

    all_iv_currs.append(iv_curr)

import pdb
pdb.set_trace()

axs[0].plot(dat['engine.time'], dat['membrane.V'], 'k')

for ax in axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axs[0].set_ylabel('Voltage (mV)')
axs[1].set_ylabel(r'$I_m$ (pA/pF)')
axs[1].set_xlabel('Time (ms)')

axs[1].legend()

plt.show()

fig, ax = plt.subplots(1, 1, figsize=(12,8))

for i, iv_currs in enumerate(all_iv_currs):
    ax.plot(list(range(-120, 60, 10)), iv_currs, label=labs[i], marker='o')


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_ylabel('Current (pA/pF)')
ax.set_xlabel('Voltage (mV)')
plt.legend()
plt.show()
