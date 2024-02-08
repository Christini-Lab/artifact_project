from deap import base, creator, tools
import myokit
import random
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
import pickle
from utility_classes import VCProtocol, VCSegment
from os import mkdir, listdir

# 0. Add expeirmental artifact to Kernik-Clancy
# 1. Create an object for voltage clamp protocols
# 2. Use DEAP to create a GA
    #Initialize
    #

random.seed(10)
global C_M 
C_M = 60

global GLOBAL_SIM
mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')

p = mod.get('engine.pace')
p.set_binding(None)

c_m = mod.get('artifact.c_m')
c_m.set_rhs(C_M)

v_cmd = mod.get('artifact.v_cmd')
v_cmd.set_rhs(0)
v_cmd.set_binding('pace') # Bind to the pacing mechanism
holding_proto = myokit.Protocol()
holding_proto.add_step(-80, 30000)
t = holding_proto.characteristic_time()
GLOBAL_SIM = myokit.Simulation(mod, holding_proto)
dat = GLOBAL_SIM.run(t)


def run_ga(toolbox):
    """
    Runs an instance of the genetic algorithm.

    Returns
    -------
        final_population : List[Individuals]
    """
    print('Evaluating initial population.')

    # 3. Calls _initialize_individuals GA_CONFIG.population size number of times and returns the initial population
    population = toolbox.population(GA_CONFIG.population_size)
    #population[0][0].plot_protocol(is_shown=True)

    # 4. Calls _evaluate_fitness on every individual in the population
    fitnesses = toolbox.map(toolbox.evaluate, population)
    for ind, fit in zip(population, fitnesses):
        ind.fitness.values = (fit,)

    # Note: visualize individual fitnesses with: population[0].fitness
    gen_fitnesses = [ind.fitness.values[0] for ind in population]

    print(f'\tAvg fitness is: {np.mean(gen_fitnesses)}')
    print(f'\tBest fitness is {np.max(gen_fitnesses)}')

    # Store initial population details for result processing.
    final_population = [population]

    for generation in range(1, GA_CONFIG.max_generations):
        print('Generation {}'.format(generation))
        # Offspring are chosen through tournament selection. They are then
        # cloned, because they will be modified in-place later on.

        # 5. DEAP selects the individuals 
        selected_offspring = toolbox.select(population, len(population))

        offspring = [toolbox.clone(i) for i in selected_offspring]

        # 6. Mate the individualse by calling _mate()
        for i_one, i_two in zip(offspring[::2], offspring[1::2]):
            if random.random() < GA_CONFIG.mate_probability:
                toolbox.mate(i_one, i_two)
                del i_one.fitness.values
                del i_two.fitness.values

        # 7. Mutate the individualse by calling _mutate()
        for i in offspring:
            if random.random() < GA_CONFIG.mutate_probability:
                toolbox.mutate(i)
                del i.fitness.values

        # All individuals who were updated, either through crossover or
        # mutation, will be re-evaluated.
        
        # 8. Evaluating the offspring of the current generation
        updated_individuals = [i for i in offspring if not i.fitness.values]
        fitnesses = toolbox.map(toolbox.evaluate, updated_individuals)
        for ind, fit in zip(updated_individuals, fitnesses):
            ind.fitness.values = (fit,)

        population = offspring

        gen_fitnesses = [ind.fitness.values[0] for ind in population]

        print(f'\tAvg fitness is: {np.mean(gen_fitnesses)}')
        print(f'\tBest fitness is {np.max(gen_fitnesses)}')

        final_population.append(population)

    return final_population


def _initialize_individual():
    # Builds a list of parameters using random upper and lower bounds.
    segments = []

    for i in range(0, GA_CONFIG.num_segments):
        if GA_CONFIG.steps_only:
            is_ramp = 0
        else:
            is_ramp = random.randint(0, 1)
        if  is_ramp == 0:
            #step
            duration = random.uniform(GA_CONFIG.vc_min_duration,
                    GA_CONFIG.vc_max_duration)
            step_voltage = random.uniform(GA_CONFIG.vc_min_voltage,
                    GA_CONFIG.vc_max_voltage)
            segment = VCSegment(duration=duration, start_voltage=step_voltage)
        else:
            #ramp
            duration = random.uniform(GA_CONFIG.vc_min_duration,
                    GA_CONFIG.vc_max_duration)
            start_voltage = random.uniform(GA_CONFIG.vc_min_voltage,
                    GA_CONFIG.vc_max_voltage)
            end_voltage = random.uniform(GA_CONFIG.vc_min_voltage,
                    GA_CONFIG.vc_max_voltage)
            segment = VCSegment(duration=duration,
                                start_voltage=start_voltage,
                                end_voltage=end_voltage)

        segments.append(segment)

    proto = VCProtocol(segments)

    return proto


def _evaluate_fitness(ind):
    # get simulation results
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
            as_array = np.abs(np.array(dat[curr])) / C_M
            tot_current += as_array
        else:
            as_array = np.abs(np.array(dat[curr]))
            tot_current += as_array

    if curr == 'artifact.i_leak':
        contrib = np.abs(np.array(dat[GA_CONFIG.target_current]) /
                                                            C_M) / tot_current
    else:
        contrib = np.abs(np.array(dat[GA_CONFIG.target_current])) / tot_current
    if GA_CONFIG.avg_window > 1:
        w = GA_CONFIG.avg_window
        max_isolation = np.max(np.convolve(contrib, np.ones(w), 'valid') / w)
    else:
        max_isolation = np.max(contrib)
    #max_idx = np.argmax(np.abs(np.array(dat[GA_CONFIG.target_current])) / tot_current)
    #fig, axs = plt.subplots(4, 1, sharex=True, figsize=(12, 8))
    #axs[0].plot(dat['engine.time'], dat['membrane.V'])
    #axs[1].plot(dat['engine.time'], np.array(dat['membrane.i_ion']) / C_M)
    #axs[2].plot(dat['engine.time'], np.array(dat['artifact.i_leak']) / C_M)
    #axs[3].plot(dat['engine.time'], tot_current)
    #plt.show()
    return max_isolation


def simulate_model(ind):
    mod = myokit.load_model('mmt-files/kernik_2019_NaL_art.mmt')

    p = mod.get('engine.pace')
    p.set_binding(None)

    c_m = mod.get('artifact.c_m')
    c_m.set_rhs(C_M)

    v_cmd = mod.get('artifact.v_cmd')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-80, 30000)
    mod.set_state(GLOBAL_SIM.state())

    # Get protocol to run
    piecewise_function, segment_dict, t_max = ind[0].get_myokit_protocol()
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


class GAConfiguration():
    def __init__(self,
                 population_size,
                 max_generations,
                 vc_max_duration,
                 vc_min_duration,
                 vc_max_voltage,
                 vc_min_voltage,
                 num_segments,
                 target_current,
                 mate_probability,
                 mutate_probability,
                 tournament_size,
                 steps_only,
                 avg_window
                 ):
        self.population_size = population_size
        self.max_generations = max_generations
        self.vc_max_duration = vc_max_duration
        self.vc_min_duration = vc_min_duration
        self.vc_max_voltage = vc_max_voltage
        self.vc_min_voltage = vc_min_voltage
        self.num_segments = num_segments
        self.target_current = target_current
        self.mate_probability = mate_probability
        self.mutate_probability = mutate_probability
        self.tournament_size = tournament_size
        self.steps_only = steps_only
        self.avg_window = avg_window

    def write_meta(self, path):
        f = open(f"{path}/meta.txt","w+")
        f.write(f"Population size: {self.population_size}\n")
        f.write(f"Max generations: {self.max_generations}\n")
        f.write(f"Window size: {self.avg_window}\n")
        f.write(f"Steps only: {self.steps_only}\n")
        f.write(f"Number of segments: {self.num_segments}\n")
        f.write(f"Mate probability: {self.mate_probability}\n")
        f.write(f"Mutate probability: {self.mutate_probability}\n")
        f.write(f"Tournament size: {self.tournament_size}\n")
        f.write(f"Duration: ({self.vc_min_duration}, {self.vc_max_duration})\n")
        f.write(f"Voltage: ({self.vc_min_voltage}, {self.vc_max_voltage})\n")
        f.close()


def _mate(i_one, i_two):
    i_one = i_one[0]
    i_two = i_two[0]

    rand_steps = [*range(0, len(i_one.segments))]
    random.shuffle(rand_steps)

    for i in range(len(i_one.segments)):
        if random.random() < GA_CONFIG.mate_probability:
            i_one.segments[i], i_two.segments[rand_steps[i]] = (i_two.segments[rand_steps[i]], i_one.segments[i])


def _mutate(ind):
    ind = ind[0]
    for seg in ind.segments:
        num = random.random()

        if num < GA_CONFIG.mutate_probability:
            is_not_accepted_mutation = True 
            j = 0
            while is_not_accepted_mutation:
                j += 1
                dur = random.gauss(seg.duration, 50)
                start_v = random.gauss(seg.start_voltage, 5)
                if seg.end_voltage is not None:
                    end_v = random.gauss(seg.end_voltage, 5)
                else:
                    end_v = 0

                if ((dur < GA_CONFIG.vc_min_duration) or
                    (dur > GA_CONFIG.vc_max_duration) or 
                    (start_v < GA_CONFIG.vc_min_voltage) or 
                    (end_v > GA_CONFIG.vc_max_voltage)):
                    is_not_accepted_mutation = True
                else:
                    is_not_accepted_mutation = False
                    seg.duration = dur
                    seg.start_voltage = start_v
                    if seg.end_voltage is not None:
                        seg.end_voltage = end_v


def start_ga(target_current):
    # 1. Initializing GA hyperparameters
    creator.create('FitnessMax', base.Fitness, weights=(1.0,))

    creator.create('Individual', list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    toolbox.register('init_param',
                     _initialize_individual)
    toolbox.register('individual',
                     tools.initRepeat,
                     creator.Individual,
                     toolbox.init_param,
                     n=1)
    toolbox.register('population',
                     tools.initRepeat,
                     list,
                     toolbox.individual)

    toolbox.register('evaluate', _evaluate_fitness)
    toolbox.register('select',
                     tools.selTournament,
                     tournsize=GA_CONFIG.tournament_size)
    toolbox.register('mate', _mate)
    toolbox.register('mutate', _mutate)

    # To speed things up with multi-threading
    p = Pool()
    toolbox.register("map", p.map)

    # Use this if you don't want multi-threading
    #toolbox.register("map", map)

    # 2. Calling the GA to run
    final_population = run_ga(toolbox)

    #last_gen = final_population[-1]
    #last_gen.sort(key=lambda x: x.fitness.values[0], reverse=true)
    #last_gen[0][0].plot_with_curr(ga_config.target_current)

    return final_population


def main():
    curr = 'artifact.i_leak'
    global GA_CONFIG
    GA_CONFIG = GAConfiguration(population_size=100,
                                      max_generations=10,
                                      vc_max_duration=1000,
                                      vc_min_duration=10,
                                      vc_max_voltage=60,
                                      vc_min_voltage=-120,
                                      num_segments=4,
                                      target_current=curr,
                                      mate_probability=0.2,
                                      mutate_probability=0.2,
                                      tournament_size=2, 
                                      steps_only=False,
                                      avg_window=10)

    pop = start_ga(curr)
    all_folders = listdir('results')
    all_folders = [f for f in all_folders if 'exp' in f]
    mkdir(f'./results/exp_{len(all_folders)}')
    dump_folder = f'results/exp_{len(all_folders)}'
    GA_CONFIG.write_meta(dump_folder)

    pickle.dump(pop, open(f"{dump_folder}/{curr.split('.')[0]}.pkl", "wb" ))


if __name__ == '__main__':
    main()

