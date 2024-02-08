import pickle
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from utility_classes import VCProtocol, VCSegment
from os import listdir
import myokit
from deap import base, creator, tools



class MyCustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        try:
            creator.create('FitnessMax', base.Fitness, weights=(1.0,))
            creator.create('Individual', list, fitness=creator.FitnessMax)
            if module == '__main__':
                module = 'vc_opt_ga'
            return super().find_class(module, name)
        except:
            import pdb
            pdb.set_trace()


def plot_with_curr_contribution(ind, target_curr, is_shown=False, cm=60):
    cm=60
    dat = simulate_model(ind)

    if dat is None:
        print("Error in model run")
        return 0

    curr_names = ['ik1.i_K1',
                  'ito.i_to',
                  'ikr.i_Kr',
                  'iks.i_Ks',
                  'ical.i_CaL',
                  'icat.i_CaT',
                  'inak.i_NaK',
                  'ina.i_Na',
                  'inaca.i_NaCa',
                  'ipca.i_PCa',
                  'ifunny.i_f',
                  'ibna.i_b_Na',
                  'ibca.i_b_Ca',
                  'INaL.INaL',
                  'artifact.i_leak']

    tot_current = np.zeros(len(dat[curr_names[0]]))
    for curr in curr_names:
        if curr == 'artifact.i_leak':
            as_array = np.abs(np.array(dat[curr])) / cm
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array


    contributions = np.abs(np.array(dat[target_curr])) / tot_current

    max_contrib = np.max(contributions)
    max_arg = np.argmax(contributions)

    fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))
    fig.suptitle(f'{round(max_contrib*100, 2)}% at {max_arg/10} ms', fontsize=16)
    piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
    #t_max = 2000
    times = np.arange(0, t_max, 0.1)

    axs[0].plot(times, dat['artifact.v_cmd'], 'k--')
    axs[0].plot(times, dat['membrane.V'])
    axs[1].plot(times, np.array(dat['artifact.i_out']) / cm)
    axs[2].plot(times, np.array(dat[target_curr]))
    axs[3].plot(times, contributions)

    axs[0].set_ylabel('Vcmd')
    axs[1].set_ylabel('I_out')
    axs[2].set_ylabel(target_curr)
    axs[3].set_ylabel('Contribution')
    axs[3].set_ylim(0, 1)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    if is_shown:
        plt.show()
    else:
        return fig


def plot_vc_protocols(folder='exp_2'):
    files = ['INaL.pkl', 'ikr.pkl', 'ical.pkl', 'ina.pkl', 'ito.pkl', 'ik1.pkl', 'ifunny.pkl', 'iks.pkl']
    currents = ['INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ina.i_Na', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks'] 

    for f_name, target_curr in dict(zip(files, currents)).items():
        if f_name not in listdir(f'results/{folder}'):
            continue
        ind = get_best_proto(folder, f_name)

        plot_with_curr_contribution(ind, target_curr=target_curr)

        plt.savefig(f'results/{folder}/{f_name.split(".")[0]}.pdf')


def compare_contributions(folder):
    files = ['ina.pkl', 'INaL.pkl', 'ikr.pkl', 'ical.pkl',                             'ito.pkl', 'ik1.pkl', 'ifunny.pkl', 'iks.pkl']
    currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks'] 

    long_proto = pickle.load(open(f'./results/{folder}/long_proto.pkl', 'rb'))

    for f_name, target_curr in dict(zip(files, currents)).items():
        ind = get_best_proto(folder, f_name)


        max_isolation, max_idx = get_current_contribution(ind[0], target_curr)
        long_max, long_idx = get_current_contribution(long_proto, target_curr)

        print(
          f'{target_curr}: Original max: {max_isolation} – long max: {long_max} at {long_idx/10}')
    

def plot_proto(ind, is_shown=True):
    cm=60
    dat = simulate_model(ind)

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
    times = np.arange(0, t_max, 0.1)

    axs[0].plot(times, dat['artifact.v_cmd'], 'k--')
    axs[0].plot(times, dat['membrane.V'])
    axs[1].plot(times, np.array(dat['artifact.i_out']) / cm)

    fs = 16
    axs[0].set_ylabel('mV', fontsize=fs)
    axs[1].set_ylabel('I_out (A/F)', fontsize=fs)
    axs[1].set_xlabel('Time (ms)', fontsize=fs)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    if is_shown:
        plt.show()
    else:
        return fig


def get_current_contribution(proto, target_current, window=10):
    dat = simulate_model([proto])
    capacitance = 60

    curr_names = ['ik1.i_K1',
                  'ito.i_to',
                  'ikr.i_Kr',
                  'iks.i_Ks',
                  'ical.i_CaL',
                  'icat.i_CaT',
                  'inak.i_NaK',
                  'ina.i_Na',
                  'inaca.i_NaCa',
                  'ipca.i_PCa',
                  'ifunny.i_f',
                  'ibna.i_b_Na',
                  'ibca.i_b_Ca',
                  'INaL.INaL',
                  'artifact.i_leak']

    tot_current = np.zeros(len(dat[curr_names[0]]))
    for curr in curr_names:
        if curr == 'artifact.i_leak':
            as_array = np.abs(np.array(dat[curr])) / capacitance
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array


    contrib = np.abs(np.array(dat[target_current])) / tot_current

    if window > 1:
        smoothed_contrib = (np.convolve(contrib,
                                np.ones(window), 'valid') / window)
        max_isolation = np.max(smoothed_contrib)
        max_idx = np.argmax(smoothed_contrib)

    else:
        max_isolation = np.max(contrib)
        max_idx = np.argmax(contrib)

    return max_isolation, max_idx 


def simulate_model(proto):
    mod = myokit.load_model('./mmt-files/kernik_2019_NaL_art.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    c_m = mod.get('artifact.c_m')
    c_m.set_rhs(60)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)
    t = holding_proto.characteristic_time()
    holding_sim = myokit.Simulation(mod, holding_proto)
    dat = holding_sim.run(t)
    mod.set_state(holding_sim.state())

    # Get protocol to run
    piecewise_function, segment_dict, t_max = proto.get_myokit_protocol()
    mem = mod.get('artifact')

    for v_name, st in segment_dict.items():
        v_new = mem.add_variable(v_name)
        v_new.set_rhs(st)

    vp = mem.add_variable('vp')
    vp.set_rhs(0)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_binding(None)
    vp.set_binding('pace')

    v_cmd.set_rhs(piecewise_function)
    times = np.arange(0, t_max, 0.1)
    #times = np.arange(0, 2000, 0.1)
    ## CHANGE THIS FROM holding_proto TO SOMETHING ELSE
    sim = myokit.Simulation(mod, holding_proto)
    try:
        dat = sim.run(t_max, log_times=times)
    except:
        return None

    return dat


def plot_all_curr_contributions(proto, target_curr=None, is_shown=True):
    cm=60
    dat = simulate_model(proto)

    if dat is None:
        print("Error in model run")
        return 0
    
    color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'brown', 'lightsalmon',
                    'goldenrod', 'greenyellow', 'aquamarine', 'turquoise',
                    'lightsteelblue', 'rebeccapurple', 'fuchsia'] 

    curr_names = ['ik1.i_K1',
                  'ito.i_to',
                  'ikr.i_Kr',
                  'iks.i_Ks',
                  'ical.i_CaL',
                  'icat.i_CaT',
                  'inak.i_NaK',
                  'ina.i_Na',
                  'inaca.i_NaCa',
                  'ipca.i_PCa',
                  'ifunny.i_f',
                  'ibna.i_b_Na',
                  'ibca.i_b_Ca',
                  'INaL.INaL',
                  'artifact.i_leak']

    color_key = dict(zip(curr_names, color_list))

    tot_current = np.zeros(len(dat[curr_names[0]]))

    if target_curr is None:
        t_range = [0, len(dat[curr_names[0]])]
    else:
        max_cont, max_idx= get_current_contribution(ind[0], target_curr)

        t_range = [max_idx- 1000, max_idx + 1000]
    
    for curr in curr_names:
        if curr == 'artifact.i_leak':
            as_array = np.abs(np.array(dat[curr])) / cm
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array


    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    piecewise_function, segment_dict, t_max = proto.get_myokit_protocol()
    #t_max = 2000
    times = np.arange(0, t_max, 0.1)

    times = times[t_range[0]:t_range[1]]
    axs[0].plot(times, dat['artifact.v_cmd'][t_range[0]:t_range[1]], 'k--')
    axs[0].plot(times, dat['membrane.V'][t_range[0]:t_range[1]])
    axs[1].plot(times, np.array(
        dat['artifact.i_out'][t_range[0]:t_range[1]]) / cm)

    currents = ['ina.i_Na', 'INaL.INaL', 'ikr.i_Kr', 'ical.i_CaL', 'ito.i_to',
                            'ik1.i_K1', 'ifunny.i_f', 'iks.i_Ks',
                            'artifact.i_leak'] 

    for curr in curr_names:
        if curr == 'artifact.i_leak':
            contributions = np.abs(
                    np.array(dat[curr][t_range[0]:t_range[1]]) / cm) / (
                            tot_current[t_range[0]:t_range[1]] )
        else:
            contributions = np.abs(
                    np.array(dat[curr][t_range[0]:t_range[1]])) / (
                            tot_current[t_range[0]:t_range[1]])

        if np.max(contributions) > .1:
            axs[2].plot(times, contributions, label=curr, c=color_key[curr])
        else:
            axs[2].plot(times, contributions, 'k')

    if target_curr is not None:
        times = np.arange(0, t_max, 0.1)
        axs[2].axvspan(
                times[max_idx-25], times[max_idx+25],facecolor='g', alpha=0.25)


    axs[0].set_ylabel('Vcmd')
    axs[1].set_ylabel('I_out')
    axs[2].set_ylabel('Contributions')
    axs[2].set_ylim(0, 1)
    axs[2].legend()

    if target_curr is not None:
        fig.suptitle(
                f'Max contribution for {target_curr} is {round(max_cont*100, 2)}%', fontsize=18)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    if is_shown:
        plt.show()
    else:
        return fig


def test_manual_protos():
    proto = VCProtocol([VCSegment(2000, 5), VCSegment(800, -20, -100), VCSegment(2000, -60)])
    #proto = VCProtocol([VCSegment(1000, -60),
    #                    VCSegment(1000, -65),
    #                    VCSegment(1000, -70),
    #                    VCSegment(1000, -75),
    #                    VCSegment(1000, -80),
    #                    VCSegment(1000, -85),
    #                    VCSegment(1000, -90),
    #                    VCSegment(1000, -95),
    #                    VCSegment(1000, -100),
    #                    VCSegment(1000, -105)])

    #proto = VCProtocol([VCSegment(1000, -60),
    #                    VCSegment(2000, -100),
    #                    VCSegment(1000, -70)])
    #proto = VCProtocol([VCSegment(1000, -60),
    #                    VCSegment(2000, 50),
    #                    VCSegment(1000, -90)])

    plot_all_curr_contributions(proto)

    #dat = simulate_model(proto)
    #fig, axs = plt.subplots(5, 1, sharex=True, figsize=(12, 8))
    #piecewise_function, segment_dict, t_max = proto.get_myokit_protocol()
    #times = np.arange(0, t_max, 0.1)
    #cm=60

    #axs[0].plot(times, dat['artifact.v_cmd'], 'k--')
    #axs[0].plot(times, dat['membrane.V'])
    #axs[1].plot(times, np.array(dat['artifact.i_out']) / cm)
    #axs[2].plot(times, np.array(dat['ik1.i_K1']))
    #axs[3].plot(times, np.array(dat['ifunny.i_f']))
    #axs[4].plot(times, np.array(dat['artifact.i_leak']) / cm)

    #axs[1].set_ylabel('I_out')
    #axs[2].set_ylabel('I_K1')
    #axs[3].set_ylabel('I_f')
    #axs[4].set_ylabel('I_leak')

    #plt.show()


def get_best_proto(folder, f_name):
    with open(f'results/{folder}/{f_name}', 'rb') as f:
        unpickler = MyCustomUnpickler(f)
        obj = unpickler.load()

    pop = pickle.load(open(f'results/{folder}/{f_name}', "rb" ))
    last_gen = pop[-1]
    last_gen.sort(key=lambda x: x.fitness.values[0], reverse=True)

    return last_gen[0]


def run_ga_artifact_proto():
    ind = get_best_proto('/exp_0', 'artifact.pkl') 
    proto = ind[0]
    #proto.segments = proto.segments[2:]
    plot_all_curr_contributions(proto)

#run_ga_artifact_proto()
test_manual_protos()
