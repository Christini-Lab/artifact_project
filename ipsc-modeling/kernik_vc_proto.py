import myokit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm



def plot_proto_response(proto, current, with_color=True, ramps=None):
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    if ramps is not None:
        c = mod.get('membrane')

        v = c.get('V')
        v.set_binding(None)

        piecewise = 'piecewise('

        for i, r in enumerate(ramps):
            # Add a v1 variable
            v1 = c.add_variable(f'v{i}')
            v1.set_rhs(r[0])
            #v1.set_rhs('-150 + 0.1 * engine.time')

            piecewise = f'{piecewise} {r[1]}, v{i}, '

        # Add a p variable
        vp = c.add_variable('vp')
        vp.set_rhs(0)
        vp.set_binding('pace')

        # Set a new right-hand side equation for V
        v.set_rhs(f'{piecewise} vp)')

    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()

    times = np.arange(0, t, 1)

    dat = sim.run(t, log_times=times)


    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    currents = ['ik1.i_K1', 'ito.i_to', 'ikr.i_Kr', 'iks.i_Ks', 'ical.i_CaL',
                    'icat.i_CaT', 'inak.i_NaK', 'ina.i_Na', 'inaca.i_NaCa',
                    'ipca.i_PCa', 'ifunny.i_f', 'ibna.i_b_Na', 'ibca.i_b_Ca']

    contributions = []
    for idx, t in enumerate(times):
        all_currs = [np.abs(dat[c][idx]) for c in currents]
        contributions.append(np.abs(dat[current][idx])/sum(all_currs))

    axs[0].plot(dat['engine.time'], dat['membrane.V'], 'k')
    if with_color:
        axs[1].scatter(dat['engine.time'], dat['membrane.i_ion'], c=contributions,
                cmap=cm.copper)
    else:
        axs[1].scatter(dat['engine.time'], dat['membrane.i_ion'], c='k')

    fs = 22
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=fs-4)



    print(f'Max isolation of {max(contributions[4:])}%')
    plt.savefig('./figs/ical_cont.png')
    plt.show()


def plot_drug_response(proto, with_drug=True):
    dats = []
    for dr_val in [1, .2]:
        mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')

        mod['ical']['g_scale'].set_rhs(dr_val)

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
        dats.append(dat)


    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))


    axs[0].plot(dat['engine.time'], dat['membrane.V'], 'k')

    st = ['k', 'r--']
    labs = ['Baseline', r'$I_{CaL}$ Block']
    for i, d in enumerate(dats):
        axs[1].plot(d['engine.time'], d['membrane.i_ion'], st[i], label=labs[i])
        if not with_drug:
            break
        

    fs = 22
    axs[0].set_ylabel('Voltage (mV)', fontsize=fs)
    axs[1].set_ylabel('Current (pA/pF)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=fs-4)


    plt.legend(fontsize=fs)
    if with_drug:
        plt.savefig('./figs/ical_iso_drug.png')
    else:
        plt.savefig('./figs/ical_iso_no_drug.png')

    plt.show()


proto = myokit.Protocol()
ramps = None

proto.add_step(-80, 200)
proto.add_step(-40, 200)

#proto 1
#proto.add_step(-100, 200)
#proto.add_step(25, 100)
#proto.add_step(-120, 200)
#proto.add_step(-95, 200)


#proto 2
#proto.add_step(-80, 200)
#proto.add_step(-30, 200)
#proto.add_step(-120, 200)
#proto.add_step(-70, 200)
#
#ramps = [['50 - .275 * engine.time', '(engine.time >= 200 and engine.time < 400)']]


##proto 3
#proto.add_step(-115, 200)
#proto.add_step(-120, 200)
#proto.add_step(-110, 200)
#proto.add_step(-105, 200)




#
#v1.set_rhs('-150 + 0.1 * engine.time')
#(engine.time >= 300 and engine.time < 700)

#ik1.i_K1
#ito.i_to
#ikr.i_Kr
#iks.i_Ks
#ical.i_CaL
#icat.i_CaT
#inak.i_NaK
#ina.i_Na
#inaca.i_NaCa
#ipca.i_PCa
#ifunny.i_f
#ibna.i_b_Na
#ibca.i_b_Ca

#plot_proto_response(proto, current='ical.i_CaL', with_color=True, ramps=ramps)
plot_drug_response(proto, True)










