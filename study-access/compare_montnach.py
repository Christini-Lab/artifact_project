import matplotlib.pyplot as plt
import numpy as np
import myokit
from plot_all_models import get_exp_sodium_proto


#UTILITY FUNCTIONS
def get_mont_sodium_proto(scale=1000):
    proto = myokit.Protocol()
    proto.add_step(-.1*scale, 2*scale)

    for i in range(-80, 66, 5):
        if i == 0:
            proto.add_step(.1/1000*scale, .050*scale)
        else:
            proto.add_step(i/1000*scale, .050*scale)

        proto.add_step(-.1*scale, 1.95*scale)

    return proto


def get_mont_data(fi):
    with open(fi) as f:
        lines = f.readlines()

    rows_as_nums = []
    for row in lines:
        rows_as_nums.append([float(n.strip()) for n in row.split('\t') if n != '\n'])

    rows_as_nums = np.array(rows_as_nums)


    num_cols = rows_as_nums.shape[1]
    
    return rows_as_nums


def get_iv_dat(voltages, current):
    iv_dat = {}

    step_idxs = np.where(np.diff(voltages) > .005)[0]

    v_steps = voltages[step_idxs + 10]
    iv_dat['Voltage'] = v_steps

    currs = []
    for idx in step_idxs:
        #start from +3
        currs.append(np.min(current[(idx+3):(idx+23)]))

    iv_dat['Current'] = currs

    return iv_dat


def simulate_model(mod, proto, with_hold=True, sample_freq=0.0001):
    if mod.time_unit().multiplier() == .001:
        scale = 1000
    else:
        scale = 1

    p = mod.get('engine.pace')
    p.set_binding(None)

    v_cmd = mod.get('voltageclamp.Vc')
    v_cmd.set_rhs(0)
    v_cmd.set_binding('pace') # Bind to the pacing mechanism

    # Run for 20 s before running the VC protocol
    if with_hold:
        holding_proto = myokit.Protocol()
        holding_proto.add_step(-.080*scale, 30*scale)
        sim = myokit.Simulation(mod, holding_proto)
        t_max = holding_proto.characteristic_time()
        sim.run(t_max)
        mod.set_state(sim.state())

    t_max = proto.characteristic_time()
    times = np.arange(0, t_max, sample_freq*scale)
    sim = myokit.Simulation(mod, proto)

    #times = np.arange(0, t_max, .00005*scale) 
    dat = sim.run(t_max, log_times=times)

    return dat, times


def get_baseline_iv(model_name, sample_freq):
    if model_name == 'Paci':
        mod = myokit.load_model('./mmt_files/paci.mmt')
        mem_name = 'Membrane.Vm'
        ion_name = 'Membrane.i_ion'
        scale = 1
        g_na = mod['i_Na']['g_Na'].value()
        mod['i_Na']['g_Na'].set_rhs(g_na*1)
    elif model_name == 'Kernik':
        mod = myokit.load_model('./mmt_files/kernik_2019_mc.mmt')
        mem_name = 'membrane.V'
        ion_name = 'membrane.i_ion'
        scale = 1000 
        g_na = mod['ina']['g_scale'].value()
        mod['ina']['g_scale'].set_rhs(g_na*.8)
    elif model_name == 'Ord':
        mod = myokit.load_model('./mmt_files/ord_na.mmt')
        mem_name = 'Membrane.V'
        ion_name = 'Membrane.i_ion'
        scale = 1000
        g_na = mod['INa']['GNa'].value()
        mod['INa']['GNa'].set_rhs(g_na*.301)
    else:
        return

    p = mod.get('engine.pace')
    p.set_binding(None)

    v = mod.get(mem_name)
    v.demote()
    v.set_rhs(0)
    v.set_binding('pace') # Bind to the pacing mechanism

    holding_proto = myokit.Protocol()
    holding_proto.add_step(-.080*scale, 30*scale)

    t = holding_proto.characteristic_time()
    sim = myokit.Simulation(mod, holding_proto)
    dat = sim.run(t)

    mod.set_state(sim.state())

    proto = get_mont_sodium_proto(scale)
    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, sample_freq*scale)

    dat = sim.run(t, log_times=times)

    voltages = np.array(dat[mem_name])
    current = np.array(dat[ion_name]) # in nA

    if model_name == 'Paci':
        voltages *= 1000

    iv_dat = get_iv_dat(voltages, current)

    return times, dat, iv_dat



#FIGURE FUNCTIONS
def plot_montnach_cpp(f):
    # Data from fig 3 of paper is:
    # 'montnach/INabad_0.600_Rs2.0_Cm20.0'
    # 'montnach/INabad_0.600_Rs5.0_Cm20.0'
    # C++ Montnach data:
    mont_dat = get_mont_data(f)
    iv_dat = get_iv_dat(mont_dat[:, 1], mont_dat[:, 3])

    #Plot Montnach C++ IV curve
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
        
    ax.plot(iv_dat['Voltage'], iv_dat['Current'], '-o')
    ax.set_xlabel('Voltage')
    ax.set_ylabel('Current')
    plt.show()
    

    fig, axs = plt.subplots(2, 1, sharex=True)

    axs[0].plot(mont_dat[:, 0], mont_dat[:, 1])
    axs[0].plot(mont_dat[:, 0], mont_dat[:, 2])
    axs[1].plot(mont_dat[:, 0], mont_dat[:, 3])
    plt.show()


def plot_multiple_montnach_cpp(GNa='0.600'):
    # Data from fig 3 of paper is:
    # 'montnach/INabad_0.600_Rs2.0_Cm20.0'
    # 'montnach/INabad_0.600_Rs5.0_Cm20.0'
    # C++ Montnach data:
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    for Rs in ['2.0', '5.0']:
        mont_dat = get_mont_data(f'montnach/INabad_{GNa}_Rs{Rs}_Cm20.0')
        iv_dat = get_iv_dat(mont_dat[:, 1], mont_dat[:, 4])

        
        ax.plot(iv_dat['Voltage'], iv_dat['Current'], '-o', label=f'Rs={Rs}M')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Voltage (mV)')
    ax.set_ylabel('Current (nA)')
    plt.legend()
    plt.show()


def plot_multiple_mmt_montnach():
    proto = get_mont_sodium_proto()
    cm = 20

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    for rs in [.002, .005]:
        mod = myokit.load_model('./mmt_files/ord_na_mont.mmt')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        
        dat, times = simulate_model(mod, proto)

        #plot vc protocol
        #fig, axs = plt.subplots(2, 1, sharex=True)
        #axs[0].plot(times, dat['voltageclamp.Vc'])
        #axs[1].plot(times, dat['voltageclamp.Iout'])
        #plt.show()

        iv_dat = get_iv_dat(np.array(dat['voltageclamp.Vc']),
                                np.array(dat['voltageclamp.Iout']))

        ax.plot(iv_dat['Voltage'], iv_dat['Current'], '-o', label=f'Rs={rs*1000}M')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Voltage (mV)')
    ax.set_ylabel('Current (pA)')
    plt.legend()
    plt.show()


def plot_mmt_vs_cpp_iv(GNa_cpp='0.600', GNa_mmt=28):
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    proto = get_mont_sodium_proto()
    cm = 20
    for rs in [.00001, .002, .005, .02]:
        mod = myokit.load_model('./mmt_files/ord_na_mont.mmt')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        iv_dat = get_iv_dat(np.array(dat['voltageclamp.Vc']),
                                np.array(dat['voltageclamp.Iout']))

        if rs < .0005:
            ax.plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        'k-o', label=f'mmt Rs<.5M')
        else:
            ax.plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=f'mmt Rs={rs*1000}M')

    for Rs in ['2.0', '5.0']:
        mont_dat = get_mont_data(f'montnach/INabad_{GNa_cpp}_Rs{Rs}_Cm20.0')
        iv_dat = get_iv_dat(mont_dat[:, 1], mont_dat[:, 4])

        
        ax.plot(iv_dat['Voltage'], iv_dat['Current'], '-*', label=f'cpp Rs={Rs}M')


    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Voltage (mV)', fontsize=16)
    ax.set_ylabel('Current (nA)', fontsize=16)
    fig.suptitle('Validating implementation of mmt ord-montanch files', fontsize=18)
    plt.legend()
    plt.show()
    

def plot_mmt_vs_cpp_iout(GNa_cpp='0.600', GNa_mmt=28):
    fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    proto = get_mont_sodium_proto()
    cm = 20
    for rs in [.002, .005]:
        mod = myokit.load_model('./mmt_files/ord_na_mont.mmt')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        if rs < .0005:
            axs[1].plot(times, np.array(dat['voltageclamp.Iout'])/1000,
                    'k-', label=f'mmt Rs<.5M')
        else:
            axs[1].plot(times, np.array(dat['voltageclamp.Iout'])/1000,
                    '-', label=f'mmt Rs={rs*1000}M')

    axs[0].plot(times, np.array(dat['voltageclamp.Vc']), '-')

    for Rs in ['2.0', '5.0']:
        mont_dat = get_mont_data(f'montnach/INabad_{GNa_cpp}_Rs{Rs}_Cm20.0')
        iv_dat = get_iv_dat(mont_dat[:, 1], mont_dat[:, 4])

        
        axs[1].plot(mont_dat[:, 0], mont_dat[:, 4], '-', label=f'cpp Rs={Rs}M')


    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[1].set_xlabel('Time (ms)', fontsize=16)
    axs[0].set_ylabel('Voltage (mV)', fontsize=16)
    axs[1].set_ylabel('Current (nA)', fontsize=16)
    fig.suptitle('Validating implementation of mmt ord-montanch files', fontsize=18)
    plt.legend()
    plt.show()


def compare_paci_ord():
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    proto = get_mont_sodium_proto(scale=1)
    cm = 20E-12
    for fi in ['paci_cardio_mont.mmt', 'paci_cardio_lei.mmt']:
        rs = 2E6
        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        #mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        voltages = np.array(dat['voltageclamp.Vc']) * 1000
        currents = np.array(dat['voltageclamp.Iout']) * 1E12

        iv_dat = get_iv_dat(voltages, currents)

        ax.plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                    '-o', label=fi)

    proto = get_mont_sodium_proto(scale=1000)
    cm = 20
    for fi in ['ord_na_mont.mmt', 'ord_na_lei.mmt']:
        rs = .002
        GNa = 22
        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        mod['INa']['GNa'].set_rhs(GNa)
        #mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        iv_dat = get_iv_dat(np.array(dat['voltageclamp.Vc']),
                                np.array(dat['voltageclamp.Iout']))

        ax.plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                    '-o', label=fi)
    plt.legend()
    plt.show()


def compare_paci_kernik_ord():
    fig, axs = plt.subplots(1, 2, figsize=(12, 8), sharey=True)

    proto = get_mont_sodium_proto(scale=1000)
    cm = 20
    for fi in ['kernik_cardio_mont.mmt', 'kernik_cardio_lei.mmt']:
        rs = .002
        GNa = .8 
        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod['geom']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        mod['ina']['g_scale'].set_rhs(GNa)
        #mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        iv_dat = get_iv_dat(np.array(dat['voltageclamp.Vc']),
                                np.array(dat['voltageclamp.Iout']))

        if 'mont' in fi:
            axs[0].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=fi)
        else:
            axs[1].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=fi)


    proto = get_mont_sodium_proto(scale=1000)
    cm = 20
    for fi in ['ord_na_mont.mmt', 'ord_na_lei.mmt']:
        rs = .002
        GNa = 22 
        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        mod['INa']['GNa'].set_rhs(GNa)
        #mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        iv_dat = get_iv_dat(np.array(dat['voltageclamp.Vc']),
                                np.array(dat['voltageclamp.Iout']))

        if 'mont' in fi:
            axs[0].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=fi)
        else:
            axs[1].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=fi)



    proto = get_mont_sodium_proto(scale=1)
    cm = 20E-12
    for fi in ['paci_cardio_mont.mmt', 'paci_cardio_lei.mmt']:
        rs = 2E6
        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod['model_parameters']['Cm'].set_rhs(cm)
        mod['voltageclamp']['cm_est'].set_rhs(cm)
        #mod['INa']['GNa'].set_rhs(GNa_mmt)
        
        dat, times = simulate_model(mod, proto)

        voltages = np.array(dat['voltageclamp.Vc']) * 1000
        currents = np.array(dat['voltageclamp.Iout']) * 1E12

        iv_dat = get_iv_dat(voltages, currents)

        if 'mont' in fi:
            axs[0].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=fi)
        else:
            axs[1].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=fi)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 16
    axs[0].set_title('Montnach Artifact', fontsize=fs+2)
    axs[1].set_title('Lei Artifact', fontsize=fs+2)
    axs[1].set_xlabel('Voltage (mV)', fontsize=fs)
    axs[0].set_xlabel('Voltage (mV)', fontsize=fs)
    axs[1].set_xlabel('Voltage (mV)', fontsize=fs)
    axs[0].set_ylabel('Current (nA)', fontsize=fs)

    plt.legend()
    plt.show()


def compare_rs_effect_on_iv(artifact_type):
    fig, axs = plt.subplots(1, 3, figsize=(16, 8), sharey=True)
    cm = 20
    comp = .7

    for i, mod_name in enumerate(['Kernik', 'Ord', 'Paci']):
        times, dat, iv_dat = get_baseline_iv(mod_name)
        axs[i].plot(iv_dat['Voltage'],
                np.array(iv_dat['Current'])*cm/1000, 'k-o', label='Ideal')

    for rseries in [.002, .005, .020]:

        for i, mod_name in enumerate(['Kernik', 'Ord', 'Paci']):
            if mod_name == 'Kernik':
                fi = f'kernik_cardio_{artifact_type}.mmt'
                scale = 1000
                rs = rseries 
                GNa = .8 
                cm_group = 'geom'
                cm_val = cm 
                gna_group = 'ina.g_scale'
            elif mod_name == 'Ord':
                fi = f'ord_na_{artifact_type}.mmt'
                scale = 1000
                rs = rseries 
                GNa = 22 
                cm_group = 'model_parameters'
                cm_val = cm
                gna_group = 'INa.GNa'
            else:
                fi = f'paci_cardio_{artifact_type}.mmt'
                scale = 1
                rs = rseries*1E9
                GNa = 1
                cm_group = 'model_parameters'
                cm_val = cm * 1E-12
                gna_group = 'i_Na.g_Na_scale'

            proto = get_mont_sodium_proto(scale=scale)
            mod = myokit.load_model(f'./mmt_files/{fi}')
            mod['voltageclamp']['rseries'].set_rhs(rs)
            if 'rseries_est' in mod['voltageclamp']._variables.keys():
                mod['voltageclamp']['rseries_est'].set_rhs(rs)
            mod[cm_group]['Cm'].set_rhs(cm_val)
            mod['voltageclamp']['cm_est'].set_rhs(cm_val)
            mod['voltageclamp']['alpha_p'].set_rhs(comp)
            mod['voltageclamp']['alpha_c'].set_rhs(comp)
            mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
            
            dat, times = simulate_model(mod, proto)

            if mod_name == 'Paci':
                voltages = np.array(dat['voltageclamp.Vc']) * 1000
                currents = np.array(dat['voltageclamp.Iout']) * 1E12
            else:
                voltages = np.array(dat['voltageclamp.Vc'])
                currents = np.array(dat['voltageclamp.Iout'])

            iv_dat = get_iv_dat(voltages, currents)

            lab = f'Rs={rseries*1E3}M'
            axs[i].plot(iv_dat['Voltage'], np.array(iv_dat['Current'])/1000,
                                        '-o', label=lab)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 16
    axs[0].set_title('Kernik', fontsize=fs+2)
    axs[1].set_title('Ord', fontsize=fs+2)
    axs[2].set_title('Paci', fontsize=fs+2)

    axs[0].set_xlabel('Voltage (mV)', fontsize=fs)
    axs[0].set_ylabel('Current (nA)', fontsize=fs)
    axs[1].set_xlabel('Voltage (mV)', fontsize=fs)
    axs[2].set_xlabel('Voltage (mV)', fontsize=fs)

    if artifact_type == 'mont':
        fig.suptitle('Montnach Artifact', fontsize=fs+4)
    else:
        fig.suptitle('Lei Artifact', fontsize=fs+4)

    print('Finished')
    plt.legend()
    plt.show()


def compare_rs_effect_on_current(mod_name, artifact_type):
    fig, axs = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    cm = 20
    comp = 0 
    sample_freq = 0.00005

    times, dat, iv_dat = get_baseline_iv(mod_name, sample_freq)
    if mod_name == 'Kernik':
        axs[0].plot(times, dat['membrane.V'], 'k-', label='Ideal')
        axs[1].plot(times, np.array(dat['membrane.i_ion'])*cm/1000,
                                                        'k-', label='Ideal')
    elif mod_name == 'Ord':
        axs[0].plot(times, dat['Membrane.V'], 'k-', label='Ideal')
        axs[1].plot(times, np.array(dat['Membrane.i_ion'])*cm/1000,
                                                    'k-', label='Ideal')
    else:
        axs[0].plot(times, dat['Membrane.Vm'], 'k-', label='Ideal')
        axs[1].plot(times, np.array(dat['Membrane.i_ion'])*cm/1000,
                                                    'k-', label='Ideal')

    for rseries in [.002, .005, .020]:
        lab = f'Rs={rseries*1000}M'
        if mod_name == 'Kernik':
            fi = f'kernik_cardio_{artifact_type}.mmt'
            scale = 1000
            rs = rseries 
            GNa = .8 
            cm_group = 'geom'
            cm_val = cm 
            gna_group = 'ina.g_scale'
        elif mod_name == 'Ord':
            fi = f'ord_na_{artifact_type}.mmt'
            scale = 1000
            rs = rseries 
            GNa = 22 
            cm_group = 'model_parameters'
            cm_val = cm
            gna_group = 'INa.GNa'
        else:
            fi = f'paci_cardio_{artifact_type}.mmt'
            scale = 1
            rs = rseries*1E9
            GNa = 1
            cm_group = 'model_parameters'
            cm_val = cm * 1E-12
            gna_group = 'i_Na.scale_Na'

        proto = get_mont_sodium_proto(scale=scale)
        mod = myokit.load_model(f'./mmt_files/{fi}')
        mod['voltageclamp']['rseries'].set_rhs(rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(rs)
        mod[cm_group]['Cm'].set_rhs(cm_val)
        mod['voltageclamp']['cm_est'].set_rhs(cm_val)
        if artifact_type == 'lei':
            mod['voltageclamp']['alpha_p'].set_rhs(comp)
            mod['voltageclamp']['alpha_c'].set_rhs(comp)
        mod[gna_group.split('.')[0]][gna_group.split('.')[1]].set_rhs(GNa)
        
        dat, times = simulate_model(mod, proto, sample_freq=sample_freq)

        if mod_name == 'Paci':
            voltages = np.array(dat['voltageclamp.Vc']) * 1000
            currents = np.array(dat['voltageclamp.Iout']) * 1E12 / 1000
        else:
            voltages = np.array(dat['voltageclamp.Vc'])
            currents = np.array(dat['voltageclamp.Iout']) / 1000

        axs[1].plot(times, currents, '-', label=lab)

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    fs = 16
    axs[0].set_ylabel('Voltage', fontsize=fs)
    axs[1].set_xlabel('Time', fontsize=fs)
    axs[1].set_ylabel('Current', fontsize=fs)

    if artifact_type == 'mont':
        fig.suptitle(f'{mod_name} – Montnach Artifact – Comp={comp*100}%',
                                                                fontsize=fs+4)
    else:
        fig.suptitle(f'{mod_name} – Lei Artifact – Comp={comp*100}%',
                                                                fontsize=fs+4)

    print('Finished')
    plt.legend()
    plt.show()


#plot_montnach_cpp('montnach/INabad_0.600_Rs5.0_Cm20.0')
#plot_multiple_montnach_cpp('2.000')
#plot_multiple_mmt_montnach()
#plot_mmt_vs_cpp_iout()
#compare_paci_ord()
#compare_paci_kernik()


#IMPORTANT
plot_mmt_vs_cpp_iv()
#compare_paci_kernik_ord()
#compare_rs_effect_on_iv('lei')
#compare_rs_effect_on_current('Paci', 'lei')
#compare_rs_effect_on_current('Paci', 'mont')
#plot_model_mult_rs('Paci', 'lei')
