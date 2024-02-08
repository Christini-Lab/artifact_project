from deap import base, creator, tools
import myokit
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool
import pickle
from utility_classes import VCProtocol, VCSegment
from os import mkdir, listdir
import xlrd
from import_rtxi import get_exp_as_df
from scipy.signal import find_peaks



# 0. Add expeirmental artifact to Kernik-Clancy
# 1. Create an object for voltage clamp protocols
# 2. Use DEAP to create a GA
    #Initialize
    #

random.seed(2)

### Core Classes
class GAConfiguration():
    def __init__(self,
                 population_size,
                 max_generations,
                 mate_probability,
                 mutate_probability,
                 tournament_size,
                 model_name,
                 cm_est,
                 ra_est,
                 params,
                 rs_comps,
                 exp_f
                 ):
        self.population_size = population_size
        self.max_generations = max_generations
        self.mate_probability = mate_probability
        self.mutate_probability = mutate_probability
        self.tournament_size = tournament_size
        self.model_name = model_name
        self.cm_est = cm_est
        self.ra_est = ra_est
        self.params = params
        self.rs_comps=rs_comps
        self.exp_f = exp_f


    def write_meta(self, path):
        f = open(f"{path}/meta.txt","w+")
        f.write(f"Experimental file: {self.exp_f}\n")
        f.write(f"Model name: {self.model_name}\n")
        f.write(f"Cell estimated Cm: {self.cm_est}\n")
        f.write(f"Cell estimated Rseries: {self.ra_est}\n")
        f.write(f"Population size: {self.population_size}\n")
        f.write(f"Max generations: {self.max_generations}\n")
        f.write(f"Mate probability: {self.mate_probability}\n")
        f.write(f"Mutate probability: {self.mutate_probability}\n")
        f.write(f"Tournament size: {self.tournament_size}\n")

        f.write(f"Adjusted params:\n")
        for p in self.params:
            f.write(f'\t{p}\n')

        f.write(f"Fitted compensation settings:\n")
        for c in self.rs_comps:
            f.write(f'\t{c}\n')

        f.close()


class Model():
    def __init__(self,
                 model_name,
                 param_values):
        self.model_name = model_name
        self.param_values = param_values
        if model_name == 'Paci':
            self.scale = 1
        else:
            self.scale = 1000

    def simulate_model(self, proto, alpha_p, alpha_c, with_hold=True):
        if self.model_name == 'Kernik':
            fi = f'kernik_cardio_lei.mmt'
            mod = myokit.load_model(f'./mmt_files/{fi}')

            mod['voltageclamp']['rseries'].set_rhs(
                                        self.param_values['rseries'] * 1E-3)
            mod['voltageclamp']['rseries_est'].set_rhs(GA_CONFIG.ra_est * 1E-3)
            mod['geom']['Cm'].set_rhs(self.param_values['cm'])
            mod['voltageclamp']['cm_est'].set_rhs(GA_CONFIG.cm_est)
            mod['ina']['g_scale'].set_rhs(self.param_values['g_Na'])
            mod['voltageclamp']['gLeak'].set_rhs(self.param_values['gLeak'])
        elif self.model_name == 'Ord':
            fi = f'ord_na_lei.mmt'
            mod = myokit.load_model(f'./mmt_files/{fi}')

            mod['voltageclamp']['rseries'].set_rhs(
                                        self.param_values['rseries'] * 1E-3)
            mod['voltageclamp']['rseries_est'].set_rhs(GA_CONFIG.ra_est * 1E-3)
            mod['model_parameters']['Cm'].set_rhs(self.param_values['cm'])
            mod['voltageclamp']['cm_est'].set_rhs(GA_CONFIG.cm_est)
            mod['INa']['g_Na_scale'].set_rhs(self.param_values['g_Na'])
            mod['voltageclamp']['gLeak'].set_rhs(self.param_values['gLeak'])
        elif self.model_name == 'Paci':
            fi = f'paci_cardio_lei.mmt'
            mod = myokit.load_model(f'./mmt_files/{fi}')
            mod['voltageclamp']['rseries'].set_rhs(
                                        self.param_values['rseries'] * 1E6)
            mod['voltageclamp']['rseries_est'].set_rhs(GA_CONFIG.ra_est * 1E6)
            mod['model_parameters']['Cm'].set_rhs(self.param_values['cm'] * 1E-12)
            mod['voltageclamp']['cm_est'].set_rhs(GA_CONFIG.cm_est * 1E-12)
            mod['i_Na']['g_Na_scale'].set_rhs(self.param_values['g_Na'])
            mod['voltageclamp']['gLeak'].set_rhs(self.param_values['gLeak']/1E9)
        else:
            print('model_name is invalid')
            return

        if 'cprs' in self.param_values.keys():
            if self.model_name == 'Paci':
                mod['voltageclamp']['cprs'].set_rhs( self.param_values['cprs'] * 1E-12)
                mod['voltageclamp']['cprs_est'].set_rhs( self.param_values['cprs_est'] * 1E-12)
            else:
                mod['voltageclamp']['cprs'].set_rhs(self.param_values['cprs'])
                mod['voltageclamp']['cprs_est'].set_rhs(self.param_values['cprs_est'])



        sample_freq = .0001

        mod['voltageclamp']['ELeak'].set_rhs(
                self.param_values['ELeak']/1000*self.scale)
        mod['voltageclamp']['alpha_p'].set_rhs(alpha_p)
        mod['voltageclamp']['alpha_c'].set_rhs(alpha_c)

        p = mod.get('engine.pace')
        p.set_binding(None)

        v_cmd = mod.get('voltageclamp.Vc')
        v_cmd.set_rhs(0)
        v_cmd.set_binding('pace') # Bind to the pacing mechanism

        # Run for 20 s before running the VC protocol
        if with_hold:
            holding_proto = myokit.Protocol()
            holding_proto.add_step(-.080*self.scale, 30*self.scale)
            sim = myokit.Simulation(mod, holding_proto)
            t_max = holding_proto.characteristic_time()
            sim.run(t_max)
            mod.set_state(sim.state())

        t_max = proto.characteristic_time()
        times = np.arange(0, t_max, sample_freq*self.scale)
        sim = myokit.Simulation(mod, proto)

        dat = sim.run(t_max, log_times=times)

        if self.model_name == 'Paci':
            voltages = np.array(dat['voltageclamp.Vc']) * 1000
            currents = np.array(dat['voltageclamp.Iout']) / (
                                                        GA_CONFIG.cm_est * 1E-12)
            times = times * 1000
        else:
            voltages = np.array(dat['voltageclamp.Vc'])
            currents = np.array(dat['voltageclamp.Iout']) / GA_CONFIG.cm_est

        return times, voltages, currents



### UTILITY FUNCTIONS
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

    time = temp_dat['Time (s)'].values
    voltages = temp_dat['Voltage (V)'].values

    return exp_dat, time, voltages, cm, ra, rm


def get_param_bounds(params):
    param_bounds = []
    for param in params:
        if param == 'rseries':
            param_bounds.append([GA_CONFIG.ra_est * .7,
                                        GA_CONFIG.ra_est * 1.3, None])
        elif param == 'cm':
            param_bounds.append([GA_CONFIG.cm_est *.7, GA_CONFIG.cm_est * 1.3, None])
        elif param == 'gLeak':
            param_bounds.append([.1, 10, 'pow'])
        elif param == 'ELeak':
            param_bounds.append([-30, 30, None])
        elif param == 'g_Na':
            param_bounds.append([.1, 10, 'pow'])
        elif param == 'cprs':
            param_bounds.append([1, 10, 'pow'])
        elif param == 'cprs_est':
            param_bounds.append([1, 10, 'pow'])
        else:
            print(f'{param} is not a valid parameter')
            return

    return dict(zip(params, param_bounds))


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

    return dat, times


def get_exp_sodium_proto(scale=1000):
    proto = myokit.Protocol()
    proto.add_step(-.08*scale, .4*scale)

    for i in range(-70, 25, 10):
        if i == 0:
            proto.add_step(.1/1000*scale, .05*scale)
        else:
            proto.add_step(i/1000*scale, .05*scale)
        proto.add_step(-.08*scale, .4*scale)

    proto.add_step(.04*scale, .05*scale)
    proto.add_step(-.08*scale, .4*scale)

    proto.add_step(.06*scale, .05*scale)
    proto.add_step(-.08*scale, .1*scale)

    return proto


def save_pop_data(folder, pop):

    params = ['Fitness'] + GA_CONFIG.params

    mkdir(f'{folder}/gen_results')

    for i, gen in enumerate(pop):
        fitnesses = []
        all_param_values = [[] for p in params]

        all_vals = dict(zip(params, all_param_values))

        for ind in gen:
            all_vals['Fitness'].append(ind.fitness.values[0])

            [all_vals[p].append(v) for p, v in ind[0].param_values.items()]

        val_pd = pd.DataFrame(all_vals)

        val_pd.to_csv(f'{folder}/gen_results/gen{i}.csv', index=False)


### GA FUNCTIONS
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
        if bounds[2] is None:
            curr_val = random.uniform(bounds[0], bounds[1])
        else:
            curr_val = 10**random.uniform(
                                np.log10(bounds[0]), np.log10(bounds[1]))

        ind_params[param] = curr_val

    mod = Model(GA_CONFIG.model_name, ind_params)

    return mod


def _evaluate_fitness(ind):
    proto = get_exp_sodium_proto(scale=ind[0].scale)

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
             'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
             'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
             'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
             'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    sim_responses = {}

    voltage_step_starts = [3999, 8499, 12999, 17499, 21999, 26499,
                                30999, 35499, 39999, 44499, 48999, 53499]
    error = 0
    
    for rs_comp in GA_CONFIG.rs_comps:
        alpha_p, alpha_c = comp_dict[rs_comp]

        try:
            time, voltages, current = ind[0].simulate_model(proto, alpha_p, alpha_c)
        except:
            print("Error with simulation")
            return 1E6

        sim_responses[rs_comp] = np.array([time, voltages, current])

        exp_curr = EXP_DAT[rs_comp]

        for idx in voltage_step_starts:
            start_pt = idx+5 #0.4 ms after step
            end_pt = idx + 31 # 3 ms after step
            diff = (np.array(current[start_pt:end_pt]) -
                            np.array(exp_curr[start_pt:end_pt]))

            error += np.sum(np.abs(diff))

    return error


def _mate(i_one, i_two):
    i_one = i_one[0]
    i_two = i_two[0]

    for param, bounds in PARAM_BOUNDS.items():
        if random.random() < .5:
            i_one.param_values[param], i_two.param_values[param] = (
                            i_two.param_values[param], i_one.param_values[param])


def _mutate(ind):
    ind = ind[0]

    for param, val in ind.param_values.items():
        num = random.random()
        if num < GA_CONFIG.mutate_probability:
            is_not_accepted_mutation = True 
            while is_not_accepted_mutation:
                new_val = random.gauss(val, val*.1)
                if ((new_val > PARAM_BOUNDS[param][0]) and
                                            (new_val < PARAM_BOUNDS[param][1])):
                    ind.param_values[param] = new_val
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
    #toolbox.register("map", map)


    # Use this if you don't want multi-threading
    #toolbox.register("map", map)

    # 2. Calling the GA to run
    final_population = run_ga(toolbox)

    return final_population


def main():
    global EXP_DAT
    exp_f = '4_021921_2_alex_control'
    EXP_DAT, time, voltages, cm, ra, rm = get_exp_comp_data(exp_f)

    params = ['rseries', 'cm', 'gLeak', 'ELeak', 'g_Na', 'cprs', 'cprs_est']

    global GA_CONFIG
    GA_CONFIG = GAConfiguration(population_size=300,
                                      max_generations=20,
                                      mate_probability=0.2,
                                      mutate_probability=0.2,
                                      tournament_size=2,
                                      model_name='Paci',
                                      cm_est=cm,
                                      ra_est=ra,
                                      params=params,
                                      rs_comps=['rscomp_80', 'rspre_80', 'rspre_0'],
                                      exp_f=exp_f
                                      )

              
    global PARAM_BOUNDS
    PARAM_BOUNDS = get_param_bounds(params) 

    pop = start_ga()

    all_folders = listdir('fit_results')
    all_folders = [f for f in all_folders if 'exp' in f]
    mkdir(f'./fit_results/exp_{len(all_folders)}')
    dump_folder = f'fit_results/exp_{len(all_folders)}'
    GA_CONFIG.write_meta(dump_folder)

    save_pop_data(dump_folder, pop)


if __name__ == '__main__':
    main()
