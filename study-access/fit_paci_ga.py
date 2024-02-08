from deap import base, creator, tools
import myokit
import random
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
import pickle
from utility_classes import VCProtocol, VCSegment
from os import mkdir, listdir
import xlrd
from import_rtxi import get_exp_as_df



# 0. Add expeirmental artifact to Kernik-Clancy
# 1. Create an object for voltage clamp protocols
# 2. Use DEAP to create a GA
    #Initialize
    #

random.seed(2)

def get_exp_comp_data(f_name='4_021921_2_alex_control'):
    h5_file = f'./data/h5_files/{f_name}.h5'
    meta = f'./data/meta/{f_name}.xlsx'

    workbook = xlrd.open_workbook(meta)
    sheet = workbook.sheet_by_name('Sheet1')

    trial_names = []
    trial_nums = []

    for rownum in range(27, 37):
        trial = sheet.row_slice(rownum, 0, 4)
        trial_names.append(trial[0].value)
        trial_nums.append(int(trial[1].value))

    min_trial = min(trial_nums)

    for rownum in range(1, 10):
        if sheet.cell(rownum, 6).value == '':
            break
        
        artifact_params = [v.value for v in sheet.row_slice(rownum, 6, 10)]

        if artifact_params[0] > min_trial:
            break
        else:
            cm = artifact_params[1]
            ra = artifact_params[2]
            rm = artifact_params[3]

    exp_dat = {}

    for i, trial_num in enumerate(trial_nums):
        temp_dat = get_exp_as_df(h5_file, trial_num, cm)
        exp_dat[trial_names[i]] = temp_dat['Current (pA/pF)'].values

    return exp_dat, temp_dat['Time (s)'].values, cm, ra, rm


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
    print(f'\tBest fitness is {np.min(gen_fitnesses)}')

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
        print(f'\tBest fitness is {np.min(gen_fitnesses)}')

        final_population.append(population)

    return final_population


def _initialize_individual():
    # Builds a list of parameters using random upper and lower bounds.
    ind_params = {}

    for param, bounds in PARAM_BOUNDS.items():
        curr_val = random.uniform(np.log10(bounds[0]), np.log10(bounds[1]))
        ind_params[param] = 10**curr_val

    return ind_params


def _evaluate_fitness(ind):
    # get simulation results
    dat, times = simulate_model(ind)
    cm = ind[0]['model_parameters.Cm'] 

    i_out = [v/cm for v in dat['voltageclamp.Iout']]
    
    exp_data = EXP_DAT[0]['rscomp_80']
    exp_times = EXP_DAT[1]

    starts = []
    ends = []

    for v in range(0, 11):
        start =  4000 * (v + 1) + 500 * v + 5
        starts.append(start)
        ends.append(start + 30)

    error = 0 

    for i, s in enumerate(starts):
        error += np.abs(np.sum(exp_data[s:ends[i]] - i_out[s:ends[i]]))

    return error


def simulate_model(ind):
    mod = myokit.load_model('./mmt_files/paci_artifact.mmt')

    for param, val in ind[0].items():
        mod[param.split('.')[0]][param.split('.')[1]].set_rhs(val)

    p = mod.get('engine.pace')
    p.set_binding(None)

    #c_m = mod.get('artifact.c_m')
    #c_m.set_rhs(C_M)

    v_cmd = mod.get('voltageclamp.Vc')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    holding_proto = myokit.Protocol()
    holding_proto.add_step(-.080, 30)
    sim = myokit.Simulation(mod, holding_proto)
    t_max = holding_proto.characteristic_time()
    sim.run(t_max)
    mod.set_state(sim.state())

    times = np.arange(0, t_max, 0.0001)
    #times = np.arange(0, 2000, 0.1)
    ## CHANGE THIS FROM holding_proto TO SOMETHING ELSE
    sodium_proto = myokit.load_protocol('./mmt_files/sodium_proto_s_exp.mmt')
    t_max = sodium_proto.characteristic_time()
    times = np.arange(0, t_max, 0.0001)
    sim = myokit.Simulation(mod, sodium_proto)

    try:
        dat = sim.run(t_max, log_times=times)
    except:
        print('error in simulation')
        import pdb
        pdb.set_trace()

    return dat, times


class GAConfiguration():
    def __init__(self,
                 population_size,
                 max_generations,
                 mate_probability,
                 mutate_probability,
                 tournament_size,
                 ):
        self.population_size = population_size
        self.max_generations = max_generations
        self.mate_probability = mate_probability
        self.mutate_probability = mutate_probability
        self.tournament_size = tournament_size

    def write_meta(self, path):
        f = open(f"{path}/meta.txt","w+")
        f.write(f"Population size: {self.population_size}\n")
        f.write(f"Max generations: {self.max_generations}\n")
        f.write(f"Mate probability: {self.mate_probability}\n")
        f.write(f"Mutate probability: {self.mutate_probability}\n")
        f.write(f"Tournament size: {self.tournament_size}\n")
        f.close()


def _mate(i_one, i_two):
    i_one = i_one[0]
    i_two = i_two[0]

    for param, bounds in PARAM_BOUNDS.items():
        if random.random() < GA_CONFIG.mate_probability:
            i_one[param], i_two[param] = (i_two[param], i_one[param])


def _mutate(ind):
    ind = ind[0]

    for param, val in ind.items():
        num = random.random()
        if num < GA_CONFIG.mutate_probability:
            is_not_accepted_mutation = True 
            while is_not_accepted_mutation:
                new_val = random.gauss(val, val*.1)
                if ((new_val > PARAM_BOUNDS[param][0]) and
                                            (new_val < PARAM_BOUNDS[param][1])):
                    ind[param] = new_val
                    is_not_accepted_mutation = False


def start_ga():
    # 1. Initializing GA hyperparameters
    creator.create('FitnessMin', base.Fitness, weights=(-1.0,))

    creator.create('Individual', list, fitness=creator.FitnessMin)

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
    global GA_CONFIG
    GA_CONFIG = GAConfiguration(population_size=200,
                                      max_generations=12,
                                      mate_probability=0.2,
                                      mutate_probability=0.2,
                                      tournament_size=2
                                      )

    global PARAM_BOUNDS
    PARAM_BOUNDS = {'i_Na.g_Na_scale': [.1, 10],
            'voltageclamp.rseries': [5E6, 50E6],
            'model_parameters.Cm': [5E6, 50E6],
            }


    global EXP_DAT
    EXP_DAT, cm, ra, rm = get_exp_comp_data()


    pop = start_ga()

    try:
        all_folders = listdir('results')
        all_folders = [f for f in all_folders if 'exp' in f]
        mkdir(f'./results/exp_{len(all_folders)}')
        dump_folder = f'results/exp_{len(all_folders)}'
        GA_CONFIG.write_meta(dump_folder)

        pickle.dump(pop, open(f"{dump_folder}/pop.pkl", "wb" ))
    except:
        import pdb
        pdb.set_trace()


if __name__ == '__main__':
    main()
