import myokit
import numpy as np
import matplotlib.pyplot as plt
from utility_classes import VCProtocol, VCSegment



def plot_mod_artifact(proto, mod_f, ax_v, ax_i):
    dats = []

    new_params = [None,
                  None]
                  #{'geom.Cm': 40,
                  #    'voltageclamp.rseries': .02,
                  #    'voltageclamp.cprs': 4.188502860812823,
                  #    'voltageclamp.voffset_eff': 2.9729555590147863,
                  #    'voltageclamp.gLeak': 0.23719047586907285} ]

    if mod_f == 'kernik':
        mods = ['kernik_2019_mc.mmt', 'kernik_artifact.mmt']
    if mod_f == 'paci':
        mods = ['paci_ms.mmt', 'paci_artifact_ms.mmt']


    for i, m_name in enumerate(mods):
        mod = myokit.load_model(f'./mmt_files/{m_name}')

        print(i)

        if new_params[i] is not None:
            for k, v in new_params[i].items():
                par_chil = k.split('.')
                mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)

        if i == 0:
            p = mod.get('engine.pace')
            p.set_binding(None)

            v = mod.get('membrane.V')
            v.demote()
            v.set_rhs(0)
            v.set_binding('pace') # Bind to the pacing mechanism

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()

        times = np.arange(0, t, .2)

        dat = sim.run(t, log_times=times)
        dats.append(dat)

    st = ['k', 'r--']
    labs = ['Baseline', 'With Exp Artifact']

    if mod_f == 'kernik':
        cm = mod['geom']['Cm'].value()
    if mod_f == 'paci':
        cm = mod['cell']['Cm'].value()

    for i, d in enumerate(dats):
        ax_v.plot(d['engine.time'], d['membrane.V'], st[i]) 

        if i == 0:
            ax_i.plot(d['engine.time'], d['membrane.i_ion'], st[i], label=labs[i])
        else:
            i_out = [i_out/cm for i_out in d['voltageclamp.Iout']]
            ax_i.plot(d['engine.time'], i_out, st[i], label=labs[i])


def plot_artifact_response(proto):
    dats = []
    new_params = [None,
                  None]
                  #{'geom.Cm': 40,
                  #    'voltageclamp.rseries': .02,
                  #    'voltageclamp.cprs': 4.188502860812823,
                  #    'voltageclamp.voffset_eff': 2.9729555590147863,
                  #    'voltageclamp.gLeak': 0.23719047586907285} ]

    for i, m_name in enumerate(['kernik_2019_mc.mmt', 'kernik_artifact.mmt']):
        mod = myokit.load_model(f'./mmt_files/{m_name}')

        if new_params[i] is not None:
            for k, v in new_params[i].items():
                par_chil = k.split('.')
                mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)


        if i == 0:
            p = mod.get('engine.pace')
            p.set_binding(None)

            v = mod.get('membrane.V')
            v.demote()
            v.set_rhs(0)
            v.set_binding('pace') # Bind to the pacing mechanism

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()

        times = np.arange(0, t, .2)

        dat = sim.run(t, log_times=times)
        dats.append(dat)


    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))



    st = ['k', 'r--']
    labs = ['Baseline', 'With Exp Artifact']
    for i, d in enumerate(dats):
        axs[0].plot(dat['engine.time'], dat['membrane.V'], st[i]) 
        if i == 0:
            axs[1].plot(d['engine.time'], d['membrane.i_ion'], st[i], label=labs[i])
        else:
            cm = mod['geom']['Cm'].value()
            i_out = [i_out/cm for i_out in d['voltageclamp.Iout']]
            axs[1].plot(d['engine.time'], i_out, st[i], label=labs[i])


    fs = 22
    axs[0].set_ylabel('Vm (mV)', fontsize=fs)
    axs[1].set_ylabel('Measured Current (pA/pF)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=fs-4)


    plt.legend(fontsize=fs)
    plt.savefig('./figs/with_art.png')

    plt.show()


def plot_kernik_artifact(proto, ax_v, ax_i):
    dats = []

    new_params = [None,
                  {'geom.Cm': 40,
                      'voltageclamp.rseries': .02,
                      'voltageclamp.cprs': 4.188502860812823,
                      'voltageclamp.voffset_eff': 2.9729555590147863,
                      'voltageclamp.gLeak': 0.23719047586907285} ]

    for i, m_name in enumerate(['kernik_2019_mc.mmt', 'kernik_artifact.mmt']):
        mod = myokit.load_model(f'./mmt_files/{m_name}')

        if new_params[i] is not None:
            for k, v in new_params[i].items():
                par_chil = k.split('.')
                mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)


        if i == 0:
            p = mod.get('engine.pace')
            p.set_binding(None)

            v = mod.get('membrane.V')
            v.demote()
            v.set_rhs(0)
            v.set_binding('pace') # Bind to the pacing mechanism

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()

        times = np.arange(0, t, .2)

        dat = sim.run(t, log_times=times)
        dats.append(dat)

    st = ['k', 'r--']
    labs = ['Baseline', 'With Exp Artifact']
    for i, d in enumerate(dats):
        ax_v.plot(d['engine.time'], d['membrane.V'], st[i]) 

        if i == 0:
            ax_i.plot(d['engine.time'], d['membrane.i_ion'], st[i], label=labs[i])
        else:
            cm = mod['geom']['Cm'].value()
            i_out = [i_out/cm for i_out in d['voltageclamp.Iout']]
            ax_i.plot(d['engine.time'], i_out, st[i], label=labs[i])


def plot_paci_artifact(proto, ax_v, ax_i):
    dats = []

    new_params = [None,
                  {}]

    for i, m_name in enumerate(['paci.mmt', 'paci_artifact.mmt']):
        mod = myokit.load_model(f'./mmt_files/{m_name}')

        if new_params[i] is not None:
            for k, v in new_params[i].items():
                par_chil = k.split('.')
                mod[k.split('.')[0]][k.split('.')[1]].set_rhs(v)


        if i == 0:
            p = mod.get('engine.pace')
            p.set_binding(None)

            v = mod.get('Membrane.Vm')
            v.demote()
            v.set_rhs(0)
            v.set_binding('pace') # Bind to the pacing mechanism

        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()

        times = np.arange(0, t, .0002)

        dat = sim.run(t, log_times=times)
        dats.append(dat)

    st = ['k', 'r--']
    labs = ['Baseline', 'With Exp Artifact']
    for i, d in enumerate(dats):
        t = [v*1000 for v in d['engine.time']]
        v = [v*1000 for v in d['Membrane.Vm']]
        ax_v.plot(t, v, st[i]) 

        if i == 0:
            ax_i.plot(t, d['Membrane.i_ion'], st[i], label=labs[i])
        else:
            cm = mod['model_parameters']['Cm'].value()
            i_out = [i_out/cm for i_out in d['voltageclamp.Iout']]
            ax_i.plot(t, i_out, st[i], label=labs[i])


def plot_example_art_models():
    proto = myokit.Protocol()
    ramps = None

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))

    step_to = -10

    proto.add_step(-80, 200)
    proto.add_step(step_to, 200)

    plot_mod_artifact(proto, 'kernik', axs[0], axs[1])
    plot_mod_artifact(proto, 'paci', axs[0], axs[2])

    fs = 22

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=fs-4)

    axs[0].set_ylabel('Vm (mV)', fontsize=fs)
    axs[1].set_ylabel('I_out (pA/pF)', fontsize=fs)
    axs[2].set_ylabel('I_out (pA/pF)', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    plt.legend(fontsize=fs)

    plt.show()


def plot_opt_art_models():
    proto = myokit.Protocol()
    ramps = None

    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))

    step_to = -10

    proto.add_step(-.08, .2)
    proto.add_step(step_to/1000, .2)
    plot_paci_artifact(proto, axs[0], axs[2])

    proto.add_step(-80, 200)
    proto.add_step(step_to, 200)

    plot_mod_artifact(proto, 'kernik_2019_mc.mmt', axs[0], axs[1])

    fs = 22

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=fs-4)

    axs[0].set_ylabel('Vm (mV)', fontsize=fs)
    axs[1].set_ylabel('I_out (pA/pF)', fontsize=fs)
    axs[2].set_ylabel('I_out (pA/pF)', fontsize=fs)
    axs[2].set_xlabel('Time (ms)', fontsize=fs)

    plt.legend(fontsize=fs)

    plt.show()


def test_paci_artifact():
    proto = myokit.Protocol()

    fig, axs = plt.subplots(8, 1, sharex=True, figsize=(12, 8))

    step_to = -10

    proto.add_step(-80, 1.1)

    mod = myokit.load_model(f'./mmt_files/paci_artifact_ms.mmt')
    mod_k = myokit.load_model(f'./mmt_files/kernik_artifact.mmt') 

    for mod in [mod, mod_k]:
        sim = myokit.Simulation(mod, proto)
        t = proto.characteristic_time()

        times = np.arange(0, t, .05)

        d = sim.run(t, log_times=times)

        t = [v for v in d['engine.time']]
        v = [v for v in d['membrane.V']]
        axs[0].plot(t, v) 

        axs[1].plot(t, d['membrane.i_ion'])

        #cm = mod['cell']['Cm'].value()
        #i_out = [i_out/cm for i_out in d['voltageclamp.Iout']]
        axs[2].plot(t, d['voltageclamp.Iout'])

        fs = 9
        axs[3].plot(t, d['voltageclamp.Vclamp'])
        axs[3].set_ylabel('Vclamp', fontsize=fs)
        axs[4].plot(t, d['voltageclamp.Vp'])
        axs[4].set_ylabel('Vp', fontsize=fs)
        axs[5].plot(t, d['voltageclamp.Vest'])
        #axs[5].plot(t, d['voltageclamp.dVestdt'])
        axs[5].set_ylabel('Vest', fontsize=fs)
        axs[6].plot(t, d['voltageclamp.Iin'])
        axs[6].set_ylabel('Iin', fontsize=fs)
        axs[7].plot(t, d['voltageclamp.Vc_prime'])
        axs[7].set_ylabel('Vc_prime', fontsize=fs)

        for ax in axs:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.tick_params(labelsize=fs-2)

        axs[0].set_ylabel('Vm (mV)', fontsize=fs)
        axs[1].set_ylabel('I_ion(pA/pF)', fontsize=fs)
        axs[2].set_ylabel('I_out (pA/pF)', fontsize=fs)
        axs[-1].set_xlabel('Time (ms)', fontsize=fs)

    plt.legend(fontsize=fs)

    plt.show()


#UTILITY
def return_vc_proto(scale=1):
    segments = [
            VCSegment(500, -80),
            VCSegment(757, 6),
            VCSegment(7, -41),
            VCSegment(101, 8.5),
            VCSegment(500, -80),
            VCSegment(106, -81),
            VCSegment(103, -2, -34),
            VCSegment(500, -80),
            VCSegment(183, -87),
            VCSegment(102, -52, 14),
            VCSegment(500, -80),
            VCSegment(272, 54, -107),
            VCSegment(103, 60),
            VCSegment(500, -80),
            VCSegment(52, -76, -80),
            VCSegment(103, -120),
            VCSegment(500, -80),
            VCSegment(940, -120),
            VCSegment(94, -77),
            VCSegment(8.1, -118),
            VCSegment(500, -80),
            VCSegment(729, 55),
            VCSegment(1000, 48),
            VCSegment(895, 59, 28)
            ]

    new_segments = []
    for seg in segments:
        if seg.end_voltage is None:
            new_segments.append(VCSegment(seg.duration*scale, seg.start_voltage*scale))
        else:
            new_segments.append(VCSegment(seg.duration*scale,
                                          seg.start_voltage*scale,
                                          seg.end_voltage*scale))

    return VCProtocol(new_segments)


plot_example_art_models()
#test_paci_artifact()
