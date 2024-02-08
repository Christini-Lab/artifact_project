import myokit
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


def simulate_model(mod, proto, with_hold=True):
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
    times = np.arange(0, t_max, 0.0001*scale)
    sim = myokit.Simulation(mod, proto)

    dat = sim.run(t_max, log_times=times)

    return dat, times


def get_exp_sodium_proto(scale):
    proto = myokit.Protocol()
    for i, val in enumerate(range(-7, 3, 1)):
        val = val / 100
        proto.add_step(-.08*scale, .4*scale)
        if val == 0:
            val += .001

        proto.add_step(val*scale, .05*scale)

    for i, val in enumerate(range(4, 7, 2)):
        val = val / 100
        proto.add_step(-.08*scale, .4*scale)
        proto.add_step(val*scale, .05*scale)

    return proto


def plot_mod_sodium(f_name):
    mod = myokit.load_model(f_name)
    if mod.time_unit().multiplier() == .001:
        scale = 1000
    else:
        scale = 1
    
    sodium_proto = get_exp_sodium_proto(scale)
    dat, times = simulate_model(mod, sodium_proto)
    cm = mod['voltageclamp']['cm_est'].value()
    i_out = [v/cm for v in dat['voltageclamp.Iout']]

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))

    axs[0].plot(times, dat['voltageclamp.Vc'])
    axs[0].plot(times, dat['Membrane.Vm'])
    axs[1].plot(times, i_out, label='Best fit')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_ylabel('Voltage (V)')
    axs[1].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    plt.show()


def compare_paci_mods():
    files = ['./mmt_files/paci_cardio_lei.mmt',
             './mmt_files/paci_cardio_mont.mmt',
             './mmt_files/paci_na_lei.mmt',
             './mmt_files/paci_na_mont.mmt']
    labs = ['Paci w Lei', 'Paci w Mont', 'Paci Na w Lei', 'Paci Na w Mont']
    st = ['g-', 'c-', 'g--', 'c--']

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    for i, f_name in enumerate(files):
        mod = myokit.load_model(f_name)

        rs = mod['voltageclamp']['rseries'].value()
        new_rs = rs * .1
        mod['voltageclamp']['rseries'].set_rhs(new_rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(new_rs)

        if mod.time_unit().multiplier() == .001:
            scale = 1000
        else:
            scale = 1
        
        sodium_proto = get_exp_sodium_proto(scale)
        dat, times = simulate_model(mod, sodium_proto)
        cm = mod['voltageclamp']['cm_est'].value()
        i_out = [v/cm for v in dat['voltageclamp.Iout']]

        axs[1].plot(times, i_out, st[i], label=labs[i])

    axs[0].plot(times, dat['voltageclamp.Vc'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_ylabel('Voltage (V)')
    axs[1].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    axs[1].legend()
    plt.show()


def compare_tord_mods():
    files = [ './mmt_files/tord_na_lei.mmt',
             './mmt_files/tord_na_mont.mmt']
    labs = ['Tor-Ord Na w Lei', 'Tor-ord Na w Mont']
    st = ['g-', 'c-']

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(12, 8))
    for i, f_name in enumerate(files):
        mod = myokit.load_model(f_name)

        rs = mod['voltageclamp']['rseries'].value()
        new_rs = rs * 1
        mod['voltageclamp']['rseries'].set_rhs(new_rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(new_rs)

        if mod.time_unit().multiplier() == .001:
            scale = 1000
        else:
            scale = 1
        
        sodium_proto = get_exp_sodium_proto(scale)
        dat, times = simulate_model(mod, sodium_proto)
        cm = mod['voltageclamp']['cm_est'].value()
        i_out = [v/cm for v in dat['voltageclamp.Iout']]

        axs[1].plot(times, i_out, st[i], label=labs[i])

    axs[0].plot(times, dat['voltageclamp.Vc'])

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_ylabel('Voltage (V)')
    axs[1].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    axs[1].legend()
    plt.show()


def get_iv_data(mod, dat, times):
    iv_dat = {}

    cm = mod['voltageclamp']['cm_est'].value()
    i_out = [v/cm for v in dat['voltageclamp.Iout']]
    v = np.array(dat['voltageclamp.Vc'])
    step_idxs = np.where(np.diff(v) > .005)[0]

    v_steps = v[step_idxs + 10]
    iv_dat['Voltage'] = v_steps

    sample_period = times[1]
    if mod.time_unit().multiplier() == .001:
        scale = 1000
    else:
        scale = 1

    currs = []
    #for idx in step_idxs:
    #    currs.append(np.min(i_out[(idx+3):(idx+23)]))

    for idx in step_idxs:
        temp_currs = i_out[(idx+3):(idx+103)]
        x = find_peaks(-np.array(temp_currs), distance=5, width=4)

        if len(x[0]) < 1:
            currs.append(np.min(temp_currs))
        else:
            currs.append(temp_currs[x[0][0]])


    iv_dat['Current'] = currs

    return iv_dat


def get_ideal_iv_dat(model):
    if model == 'Paci':
        mod = myokit.load_model('./mmt_files/paci.mmt')
        mem_name = 'Membrane.Vm'
        ion_name = 'Membrane.i_ion'
        scale = 1
        g_na = mod['i_Na']['g_Na'].value()
        mod['i_Na']['g_Na'].set_rhs(g_na*1)
    elif model == 'Tord':
        mod = myokit.load_model('./mmt_files/tord_na.mmt')
        mem_name = 'membrane.v'
        ion_name = 'membrane.i_ion'
        scale = 1000 
        g_na = mod['INa']['GNa'].value()
        mod['INa']['GNa'].set_rhs(g_na*1)
    elif model == 'Ord':
        mod = myokit.load_model('./mmt_files/ord_na.mmt')
        mem_name = 'Membrane.V'
        ion_name = 'Membrane.i_ion'
        scale = 1000 
        g_na = mod['INa']['GNa'].value()
        mod['INa']['GNa'].set_rhs(g_na*1)
    else:
        print('DONT RECOGNIZE MODEL NAME')
        import pdb
        pdb.set_trace()

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

    proto = get_exp_sodium_proto(scale)
    sim = myokit.Simulation(mod, proto)
    t = proto.characteristic_time()
    times = np.arange(0, t, .0001*scale)

    dat = sim.run(t, log_times=times)

    iv_dat = {}

    i_out = [v for v in dat[ion_name]]
    v = np.array(dat[mem_name])

    step_idxs = np.where(np.diff(v) > .005*scale)[0]

    v_steps = v[step_idxs + 10]
    iv_dat['Voltage'] = v_steps

    currs = []
    for idx in step_idxs:
        currs.append(np.min(i_out[(idx+1):(idx+21)]))

    iv_dat['Current'] = currs

    return dat, times, iv_dat


def plot_paci_cardio_vs_na_iv():
    files = ['./mmt_files/paci_cardio_lei.mmt',
             './mmt_files/paci_cardio_mont.mmt',
             './mmt_files/paci_na_lei.mmt',
             './mmt_files/paci_na_mont.mmt']
    labs = ['Paci w Lei', 'Paci w Mont', 'Paci Na w Lei', 'Paci Na w Mont']
    st = ['g-', 'c-', 'g--', 'c--']

    fig, axs = plt.subplots(1, 2, figsize=(12, 8))
    for i, f_name in enumerate(files):
        mod = myokit.load_model(f_name)

        rs = mod['voltageclamp']['rseries'].value()
        new_rs = rs
        mod['voltageclamp']['rseries'].set_rhs(new_rs)
        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(new_rs)

        if mod.time_unit().multiplier() == .001:
            scale = 1000
        else:
            scale = 1
        
        sodium_proto = get_exp_sodium_proto(scale)
        dat, times = simulate_model(mod, sodium_proto)
        cm = mod['voltageclamp']['cm_est'].value()
        i_out = [v/cm for v in dat['voltageclamp.Iout']]
        iv_dat = get_iv_data(mod, dat, times)

        axs[0].plot(iv_dat['Voltage'], iv_dat['Current'], f'{st[i]}o', label=labs[i])
        axs[1].plot(times, i_out, st[i], label=labs[i])

    dat, times, iv_dat = get_ideal_iv_dat('Paci')
    axs[0].plot(iv_dat['Voltage'], iv_dat['Current'], 'k-o', label='Baseline')
    axs[1].plot(times, dat['Membrane.i_ion'], 'k-', label='Baseline')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_xlabel('Voltage (V)')
    axs[0].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    #axs[1].set_xlim(2.6475, 2.6655)
    axs[0].set_title("Rs=20M")
    axs[1].set_title("Step to -20 mV")
    axs[1].legend()
    plt.show()


def plot_paci_vs_tord():
    files = ['./mmt_files/paci_cardio_lei.mmt',
             './mmt_files/paci_cardio_mont.mmt',
             './mmt_files/tord_na_lei.mmt',
             './mmt_files/tord_na_mont.mmt']
    labs = ['Paci w Lei', 'Paci w Mont', 'Tor-ord Na w Lei', 'Tor-ord Na w Mont']
    st = ['g-', 'c-', 'b-', 'r-']

    fig, axs = plt.subplots(1, 2, figsize=(12, 8))
    for i, f_name in enumerate(files):
        mod = myokit.load_model(f_name)

        rs = mod['voltageclamp']['rseries'].value()
        new_rs = rs * .1
        mod['voltageclamp']['rseries'].set_rhs(new_rs)
        # REMOVE
        if 'Paci' in labs[i]:
            g_na = mod['i_Na']['g_Na'].value()
            mod['i_Na']['g_Na'].set_rhs(g_na*5)
        else:
            g_na = mod['INa']['GNa'].value()
            mod['INa']['GNa'].set_rhs(g_na*5)

        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(new_rs)
            mod['voltageclamp']['alpha_p'].set_rhs(0)
            mod['voltageclamp']['alpha_c'].set_rhs(0)

        if mod.time_unit().multiplier() == .001:
            scale = 1000
        else:
            scale = 1
        
        sodium_proto = get_exp_sodium_proto(scale)
        dat, times = simulate_model(mod, sodium_proto)
        cm = mod['voltageclamp']['cm_est'].value()
        i_out = [v/cm for v in dat['voltageclamp.Iout']]
        iv_dat = get_iv_data(mod, dat, times)

        axs[0].plot(np.array(iv_dat['Voltage'])/scale, iv_dat['Current'], f'{st[i]}o', label=labs[i])
        axs[1].plot(times/scale, i_out, st[i], label=labs[i])

    dat, times, iv_dat = get_ideal_iv_dat('Paci')
    axs[0].plot(iv_dat['Voltage'], iv_dat['Current'], 'k-o', label='Base Paci')
    axs[1].plot(times, dat['Membrane.i_ion'], 'k-', label='Baseline Paci')

    dat, times, iv_dat = get_ideal_iv_dat('Tord')
    axs[0].plot(np.array(iv_dat['Voltage'])/1000, iv_dat['Current'], 'y-o', label='Base Tord')
    axs[1].plot(times/1000, np.array(dat['membrane.i_ion']), 'y-', label='Baseline Tord')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_xlabel('Voltage (V)')
    axs[0].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_xlim(3.549, 3.56)
    axs[0].set_title("Rs=2M")
    axs[1].set_title("Step to 0mV")
    axs[1].legend()
    plt.show()


def plot_tor_ord():
    files = ['./mmt_files/tord_na_lei.mmt',
             './mmt_files/tord_na_mont.mmt']
    labs = ['Tor-ord Na w Lei', 'Tor-ord Na w Mont']
    st = ['g-', 'c-']

    fig, axs = plt.subplots(1, 2, figsize=(12, 8))
    for i, f_name in enumerate(files):
        mod = myokit.load_model(f_name)

        rs = mod['voltageclamp']['rseries'].value()
        new_rs = rs * .1
        mod['voltageclamp']['rseries'].set_rhs(new_rs)
        # REMOVE
        g_na = mod['INa']['GNa'].value()
        mod['INa']['GNa'].set_rhs(g_na*1)

        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(new_rs)
            mod['voltageclamp']['alpha_p'].set_rhs(0)
            mod['voltageclamp']['alpha_c'].set_rhs(0)

        if mod.time_unit().multiplier() == .001:
            scale = 1000
        else:
            scale = 1
        
        sodium_proto = get_exp_sodium_proto(scale)
        dat, times = simulate_model(mod, sodium_proto)
        #cm = mod['voltageclamp']['cm_est'].value()
        cm = 1
        i_out = [v/cm for v in dat['voltageclamp.Iout']]
        iv_dat = get_iv_data(mod, dat, times)

        axs[0].plot(np.array(iv_dat['Voltage']), np.array(iv_dat['Current']), f'{st[i]}o', label=labs[i])
        axs[1].plot(times, i_out, st[i], label=labs[i])

    dat, times, iv_dat = get_ideal_iv_dat('Tord')
    axs[0].plot(np.array(iv_dat['Voltage']), np.array(iv_dat['Current']), 'y-o', label='Base Tord')
    axs[1].plot(times, np.array(dat['membrane.i_ion']), 'y-', label='Baseline Tord')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_xlabel('Voltage (V)')
    axs[0].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    #axs[1].set_xlim(3549, 3560)
    axs[0].set_title("Rs=2M")
    axs[1].set_title("Step to 0mV")
    axs[1].legend()
    plt.show()


def plot_ord():
    files = ['./mmt_files/ord_na_mont.mmt']
    labs = ['Ord Na w Mont']
    st = ['g-']

    fig, axs = plt.subplots(1, 2, figsize=(12, 8))
    for i, f_name in enumerate(files):
        mod = myokit.load_model(f_name)

        rs = mod['voltageclamp']['rseries'].value()
        new_rs = rs * .1
        mod['voltageclamp']['rseries'].set_rhs(new_rs)
        # REMOVE
        g_na = mod['INa']['GNa'].value()
        mod['INa']['GNa'].set_rhs(g_na*1)

        if 'rseries_est' in mod['voltageclamp']._variables.keys():
            mod['voltageclamp']['rseries_est'].set_rhs(new_rs)
            mod['voltageclamp']['alpha_p'].set_rhs(0)
            mod['voltageclamp']['alpha_c'].set_rhs(0)

        if mod.time_unit().multiplier() == .001:
            scale = 1000
        else:
            scale = 1
        
        sodium_proto = get_exp_sodium_proto(scale)
        dat, times = simulate_model(mod, sodium_proto)
        #cm = mod['voltageclamp']['cm_est'].value()
        cm = 60 
        i_out = [v/cm for v in dat['voltageclamp.Iout']]
        iv_dat = get_iv_data(mod, dat, times)

        axs[0].plot(np.array(iv_dat['Voltage']), np.array(iv_dat['Current']), f'{st[i]}o', label=labs[i])
        axs[1].plot(times, i_out, st[i], label=labs[i])

    dat, times, iv_dat = get_ideal_iv_dat('Ord')
    axs[0].plot(np.array(iv_dat['Voltage']), np.array(iv_dat['Current']), 'y-o', label='Base Tord')
    axs[1].plot(times, np.array(dat['Membrane.i_ion']), 'y-', label='Baseline Ord')

    for ax in axs:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    axs[0].set_xlabel('Voltage (V)')
    axs[0].set_ylabel('Current (A/F)')
    axs[1].set_xlabel('Time (s)')
    #axs[1].set_xlim(3549, 3560)
    axs[0].set_title("Rs=2M")
    axs[1].set_title("Step to 0mV")
    axs[1].legend()
    plt.show()


def plot_kernik_conc_change():


def main():
    #plot_mod_sodium('./mmt_files/tord_na_lei.mmt')
    #compare_paci_mods()
    #compare_tord_mods()
    plot_paci_cardio_vs_na_iv()
    #plot_paci_vs_tord()
    #plot_tor_ord()
    #plot_ord()


if __name__ == '__main__':
    main()
