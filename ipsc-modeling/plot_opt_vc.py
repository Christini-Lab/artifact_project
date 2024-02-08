import myokit
import numpy as np
import matplotlib.pyplot as plt
from utility_classes import VCProtocol, VCSegment



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


def get_mod_response():
    mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    #v = mod.get('membrane.V')
    #v.set_rhs(0)
    #v.set_binding('pace')

    vc_proto = return_vc_proto()

    proto = myokit.Protocol()
    proto.add_step(-80, 10000)

    piecewise, segment_dict, t_max = vc_proto.get_myokit_protocol()

    v = mod.get('membrane.V')
    v.demote()
    v.set_rhs(0)

    new_seg_dict = {}
    for k, vol in segment_dict.items():
        new_seg_dict[k] = vol 

    segment_dict = new_seg_dict

    mem = mod.get('membrane')
    vp = mem.add_variable('vp')
    vp.set_rhs(0)
    vp.set_binding('pace')

    for v_name, st in segment_dict.items():
        v_new = mem.add_variable(v_name)
        v_new.set_rhs(st)

    v.set_rhs(piecewise)

    t = proto.characteristic_time()
    sim = myokit.Simulation(mod, proto)

    times = np.arange(0, t_max, 0.1)

    dat = sim.run(t_max, log_times=times)

    return times, dat


def get_paci_response():
    mod = myokit.load_model('./mmt_files/paci.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    vc_proto = return_vc_proto(scale=1/1000)

    proto = myokit.Protocol()
    proto.add_step(-.080, 10)

    piecewise, segment_dict, t_max = vc_proto.get_myokit_protocol(scale=1/1000)

    v = mod.get('Membrane.Vm')
    v.demote()
    v.set_rhs(0)

    new_seg_dict = {}
    for k, vol in segment_dict.items():
        new_seg_dict[k] = vol 

    segment_dict = new_seg_dict

    mem = mod.get('Membrane')
    vp = mem.add_variable('vp')
    vp.set_rhs(0)
    vp.set_binding('pace')

    for v_name, st in segment_dict.items():
        v_new = mem.add_variable(v_name)
        v_new.set_rhs(st)
    
    v = mod.get('Membrane.Vm')
    v.set_rhs(piecewise)

    t = proto.characteristic_time()
    #t_max = 10
    sim = myokit.Simulation(mod, proto)

    times = np.arange(0, t_max, 0.0001)

    dat = sim.run(t_max, log_times=times)

    return times, dat


def plot_proto_response():
    kernik_times, kernik_dat = get_kernik_response()
    paci_times, paci_dat = get_paci_response()


    #currents = ['ik1.i_K1', 'ito.i_to', 'ikr.i_Kr', 'iks.i_Ks', 'ical.i_CaL',
    #                'icat.i_CaT', 'inak.i_NaK', 'ina.i_Na', 'inaca.i_NaCa',
    #                'ipca.i_PCa', 'ifunny.i_f', 'ibna.i_b_Na', 'ibca.i_b_Ca']

    #contributions = []
    #for idx, t in enumerate(times):
    #    all_currs = [np.abs(dat[c][idx]) for c in currents]
    #    contributions.append(np.abs(dat[current][idx])/sum(all_currs))

    fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12, 8))

    axs[0].plot(kernik_dat['engine.time'], kernik_dat['membrane.V'], 'k')

    #I_ion
    axs[1].plot(kernik_dat['engine.time'], kernik_dat['membrane.i_ion'], label='Kernik')
    axs[1].plot(np.array(paci_dat['engine.time'])*1000, np.array(paci_dat['Membrane.i_ion']), label='Paci')

    #INaCa
    axs[2].plot(kernik_dat['engine.time'], kernik_dat['inaca.i_NaCa'], label='Kernik')
    axs[2].plot(np.array(paci_dat['engine.time'])*1000, np.array(paci_dat['i_NaCa.i_NaCa']), label='Paci')

    #Cai
    axs[3].plot(kernik_dat['engine.time'], kernik_dat['cai.Cai'], label='Kernik')
    axs[3].plot(np.array(paci_dat['engine.time'])*1000, np.array(paci_dat['calcium_dynamics.Cai']), label='Paci')

    #Nai
    axs[4].plot(kernik_dat['engine.time'], kernik_dat['nai.Nai'], label='Kernik')
    axs[4].plot(np.array(paci_dat['engine.time'])*1000, np.array(paci_dat['sodium_dynamics.Nai']), label='Paci')

    axs[1].axhline(y=0, color='grey')
    axs[2].axhline(y=0, color='grey')


    fs = 14 
    axs[0].set_ylabel('mV', fontsize=fs)
    axs[1].set_ylabel(r'$I_m$ (pA/pF)', fontsize=fs)
    axs[2].set_ylabel(r'$I_{NaCa}$ (pA/pF)', fontsize=fs)
    axs[3].set_ylabel(r'$Ca_i$', fontsize=fs)
    axs[4].set_ylabel(r'$Na_i$', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelsize=fs-4)

    plt.legend()

    plt.show()


plot_proto_response()
