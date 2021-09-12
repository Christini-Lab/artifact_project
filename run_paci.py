import myokit
import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(12, 8))

proto = myokit.load_protocol('./mmt_files/sodium_proto_s.mmt')
mod = myokit.load_model('./mmt_files/paci.mmt')

p = mod.get('engine.pace')
p.set_binding(None)

# Get membrane potential, demote to an ordinary variable
#if is_v_demoted:
v = mod.get('Membrane.Vm')
v.demote()
v.set_rhs(0)
v.set_binding('pace') # Bind to the pacing mechanism

sim = myokit.Simulation(mod, proto)
t_max = proto.characteristic_time()
times = np.arange(0, t_max, 0.0001)

dat = sim.run(t_max, log_times=times)

axs[1].plot(dat['engine.time'], dat['Membrane.i_ion'], 'k--', label='Ideal')
axs[2].plot(dat['engine.time'], dat['i_Na.i_Na'], 'k--')


for alpha_p in [0, .2, .4, .6, .8]:
    mod = myokit.load_model('./mmt_files/paci_artifact.mmt')
    mod['voltageclamp']['alpha_c'].set_rhs(alpha_p)
    sim = myokit.Simulation(mod, proto)
    t_max = proto.characteristic_time()
    times = np.arange(0, t_max, 0.0001)

    dat = sim.run(t_max, log_times=times)


    axs[0].plot(dat['engine.time'], dat['Membrane.Vm'], label=f'AlphaP = {alpha_p}')
    i_out = [v/9E-11 for v in dat['voltageclamp.Iout']] 
    axs[1].plot(dat['engine.time'], i_out)
    axs[2].plot(dat['engine.time'], dat['i_Na.i_Na'])

axs[0].plot(dat['engine.time'], dat['voltageclamp.Vc'])
axs[0].legend()
plt.show()
