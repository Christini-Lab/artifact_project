import pickle
from cell_models import protocols


proto = pickle.load(open(f'data/shortened_trial_steps_ramps_200_50_4_-120_60_500_artefact_True_short.pkl', 'rb'))

piecewise_function = 'piecewise('
start_time = 0

all_segments = {}

for i, st in enumerate(proto.steps):
    new_seg = []
    
    if isinstance(st, protocols.VoltageClampStep):
        new_seg.append(f'{st.voltage}')
    else:
        slope = (st.voltage_end - st.voltage_start) / st.duration
        y_int = st.voltage_start - slope * start_time

        new_seg.append(f'{slope} * engine.time + {y_int}')

    end_time = start_time + st.duration
    new_seg.append(start_time)
    new_seg.append(end_time)

    start_time += st.duration

    all_segments[f'v{i}'] = new_seg

pickle.dump(all_segments, open('./data/proto.pkl', 'wb'))
