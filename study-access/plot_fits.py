import myokit
from os import mkdir, listdir
import pandas as pd
import xlrd
from import_rtxi import get_exp_as_df
from fit_models_to_exp import get_exp_sodium_proto
import numpy as np
import matplotlib.pyplot as plt
from compare_models_exp import get_modexp_curr_voltage_dat
from scipy.signal import find_peaks



#UTILITY CLASS
class Model():
    def __init__(self,
                 model_name,
                 cm_est,
                 ra_est,
                 param_values):
        self.model_name = model_name
        self.cm_est = cm_est
        self.ra_est = ra_est
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
            mod['voltageclamp']['rseries_est'].set_rhs(self.ra_est * 1E-3)
            mod['geom']['Cm'].set_rhs(self.param_values['cm'])
            mod['voltageclamp']['cm_est'].set_rhs(self.cm_est)
            mod['ina']['g_scale'].set_rhs(self.param_values['g_Na'])
            mod['voltageclamp']['gLeak'].set_rhs(self.param_values['gLeak'])
        elif self.model_name == 'Ord':
            fi = f'ord_na_lei.mmt'
            mod = myokit.load_model(f'./mmt_files/{fi}')

            mod['voltageclamp']['rseries'].set_rhs(
                                        self.param_values['rseries'] * 1E-3)
            mod['voltageclamp']['rseries_est'].set_rhs(self.ra_est * 1E-3)
            mod['model_parameters']['Cm'].set_rhs(self.param_values['cm'])
            mod['voltageclamp']['cm_est'].set_rhs(self.cm_est)
            mod['INa']['g_Na_scale'].set_rhs(self.param_values['g_Na'])
            mod['voltageclamp']['gLeak'].set_rhs(self.param_values['gLeak'])
        elif self.model_name == 'Paci':
            fi = f'paci_cardio_lei.mmt'
            mod = myokit.load_model(f'./mmt_files/{fi}')
            mod['voltageclamp']['rseries'].set_rhs(
                                        self.param_values['rseries'] * 1E6)
            mod['voltageclamp']['rseries_est'].set_rhs(self.ra_est * 1E6)
            mod['model_parameters']['Cm'].set_rhs(self.param_values['cm'] * 1E-12)
            mod['voltageclamp']['cm_est'].set_rhs(self.cm_est * 1E-12)
            mod['i_Na']['g_Na_scale'].set_rhs(self.param_values['g_Na'])
            mod['voltageclamp']['gLeak'].set_rhs(self.param_values['gLeak']/1E9)
        else:
            print('model_name is invalid')
            return

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
                                                        self.cm_est * 1E-12)
            times = times * 1000
        else:
            voltages = np.array(dat['voltageclamp.Vc'])
            currents = np.array(dat['voltageclamp.Iout']) / self.cm_est

        return times, voltages, currents


#UTILITY FUNCTIONS
def get_fit_information(exp_folder):
    all_gen_folders = listdir(f'{exp_folder}/gen_results')
    max_f = max([int(f.split('.')[0].split('gen')[1]) for f in all_gen_folders])
    final_pop = pd.read_csv(f'{exp_folder}/gen_results/gen{max_f}.csv')

    #get model name
    with open(f'{exp_folder}/meta.txt') as f:
        lines = f.readlines()

    model_name = lines[1].split(' ')[-1].split('\n')[0]
    f_name = lines[0].split(' ')[-1].split('\n')[0]

    f.close()

    #TODO: pull the following file name from the meta.txt file
    exp_dat, time, voltages, cm_est, rs_est, rm = get_exp_comp_data(f_name)

    target_settings = [v.split('\n')[0].split('\t')[1] for v in lines if
                                            (('rscomp_' in v) or 'rspre_' in v)]
    exp_targets = {}
    for k in target_settings:
        exp_targets[k] = exp_dat[k]

    return final_pop, exp_targets, model_name, cm_est, rs_est, f_name


def get_exp_comp_data(f_name):
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


def get_current_peak_traces(voltages, current):
    iv_dat = {}

    step_idxs = np.where(np.diff(voltages) > .005)[0]

    v_steps = voltages[step_idxs + 10]

    currs = []
    for idx in step_idxs:
        temp_currs = current[(idx-20):(idx+150)]

        currs.append(temp_currs)

    iv_dat = dict(zip(v_steps, currs))

    return iv_dat


def get_iv_dat(voltages, current):
    iv_dat = {}

    step_idxs = np.where(np.diff(voltages) > .005)[0]

    v_steps = voltages[step_idxs + 10]
    iv_dat['Voltage'] = v_steps

    currs = []
    for idx in step_idxs:
        #start from +3
        temp_currs = current[(idx+3):(idx+103)]
        x = find_peaks(-np.array(temp_currs), distance=5, width=3)
        #if len(x[0]) > 1:

        if len(x[0]) < 1:
            currs.append(np.min(temp_currs))
        else:
            currs.append(temp_currs[x[0][0]])

    iv_dat['Current'] = currs

    return iv_dat



#FIGURE FUNCTIONS
def plot_best_vs_exp(folder):
    final_pop, exp_targs, mod_name, cm_est, rs_est, f_name = get_fit_information(folder)

    # create optimal model
    best_ind = final_pop.iloc[final_pop['Fitness'].idxmin()]
    param_values = best_ind.to_dict()
    for k, v in param_values.items():
        print(f'{k}: {round(v, 2)}')
    param_values.pop('Fitness')

    if mod_name == 'Paci':
        proto = get_exp_sodium_proto(1)
    else:
        proto = get_exp_sodium_proto()

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    for k, exp_curr in exp_targs.items():
        alpha_p, alpha_c = comp_dict[k][0], comp_dict[k][1]

        mod_original = Model(mod_name, cm_est, rs_est,
                            {'rseries': rs_est,
                             'cm': cm_est,
                             'gLeak': 1,
                             'ELeak': 0,
                             'g_Na': 1})

        times, voltages, current = mod_original.simulate_model(proto, alpha_p, alpha_c)

        trace_original = get_current_peak_traces(voltages, current)


        best_mod = Model(mod_name, cm_est, rs_est, param_values)
        times, voltages, current = best_mod.simulate_model(
                                                    proto, alpha_p, alpha_c)

        exp_traces = get_current_peak_traces(voltages, exp_curr)

        mod_traces = get_current_peak_traces(voltages, current)

        fig, axs = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(12, 9))

        it = 0
        labs = ['Experimental', 'No fit', 'Best Fit']
        t = np.linspace(0, 16.9, 170)

        for it, traces in enumerate([exp_traces, trace_original, mod_traces]):
            i_fig = 0
            for voltage, curr_trace in traces.items():
                row = int(i_fig/3)
                col = np.mod(i_fig, 3)

                if voltage in [-70, -.07, -10, -.01, 10, .01]:
                    continue
                if it == 0:
                    axs[row, col].plot(t, curr_trace, 'k', label=labs[it])
                elif it == 1:
                    axs[row, col].plot(t, curr_trace, 'lightcoral', linestyle='--', label=labs[it])
                else:
                    axs[row, col].plot(t, curr_trace, 'lightcoral', label=labs[it])

                if it == 1:
                    axs[row, col].set_title(f'Step to: {voltage} mV')

                i_fig += 1


        for row in axs:
            row[0].set_ylabel('Current (A/F)')
            for ax in row:
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.yaxis.grid()
                ax.set_ylim(top=15)
                ax.axvspan(2.4, 5, alpha=.2)
                ax.set_xlim(0, 9)

        for col in [0, 1, 2]:
            axs[2, col].set_xlabel('Time (ms)')

        axs[2, 2].legend()
        fig.suptitle(
                f'{mod_name} fit to Exp – Pred: {int(100*alpha_p)}% and Comp: {int(100*alpha_c)}%',
                fontsize=18)
        
        exp_contents = listdir(folder)
        if 'trace_comparison' not in exp_contents:
            mkdir(f'{folder}/trace_comparison')
        plt.savefig(f'{folder}/trace_comparison/{k}')
        plt.show()


def compare_target_iv_curves(folder):
    final_pop, exp_targs, mod_name, cm_est, rs_est, f = get_fit_information(folder)

    # create optimal model
    best_ind = final_pop.iloc[final_pop['Fitness'].idxmin()]
    param_values = best_ind.to_dict()
    param_values.pop('Fitness')

    if mod_name == 'Paci':
        proto = get_exp_sodium_proto(1)
    else:
        proto = get_exp_sodium_proto()

    comp_dict = {'rspre_0': [0, 0], 'rspre_20': [.2, 0],
                 'rspre_40': [.4, 0], 'rspre_60': [.6, 0],
                 'rspre_80': [.8, 0], 'rscomp_0': [.7, 0],
                 'rscomp_20': [.7, .2], 'rscomp_40': [.7, .4],
                 'rscomp_60': [.7, .6], 'rscomp_80': [.7, .8]}

    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(12, 6))

    i = 0

    for k, exp_curr in exp_targs.items():
        trace_original, exp_dat = get_modexp_curr_voltage_dat(
                                                        f, mod_name, k)

        alpha_p, alpha_c = comp_dict[k][0], comp_dict[k][1]

        best_mod = Model(mod_name, cm_est, rs_est, param_values)
        mod_original = Model(mod_name, cm_est, rs_est, 
                                {'rseries': rs_est,
                                 'cm': cm_est,
                                 'gLeak': 1,
                                 'ELeak': 0,
                                 'g_Na': 1})

        times, voltages, current = best_mod.simulate_model(
                                                    proto, alpha_p, alpha_c)
        orig_times, orig_voltages, orig_current = mod_original.simulate_model(
                                                        proto, alpha_p, alpha_c)

        #GET IV DATA
        exp_iv = get_iv_dat(voltages, exp_curr)
        mod_iv = get_iv_dat(voltages, current)
        orig_iv = get_iv_dat(orig_voltages, orig_current)

        axs[i].plot(exp_iv['Voltage'], exp_iv['Current'], 'k-o', label='Exp')
        axs[i].plot(mod_iv['Voltage'], mod_iv['Current'], 'lightcoral', linestyle='-', marker='o', label='Best fit')
        axs[i].plot(orig_iv['Voltage'], orig_iv['Current'], 'lightcoral', linestyle='--', marker='o', label='No fit')

        axs[i].set_title(f'Pred: {int(100*alpha_p)}% & Comp: {int(100*alpha_c)}%')

        i += 1

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel('Voltage (mV)')

    axs[0].set_ylabel('Current (A/F)')

    axs[2].legend()
    fig.suptitle(
            f'{mod_name} fit to Exp',
            fontsize=18)
    plt.show()



#plot_best_vs_exp(folder='fit_results/exp_11')

compare_target_iv_curves(folder='fit_results/exp_11')
