import pickle
from cell_models import protocols


proto = pickle.load(open(f'data/shortened_trial_steps_ramps_200_50_4_-120_60_500_artefact_True_short.pkl', 'rb'))

piecewise_function = 'piecewise('
start_time = 0

all_segments = {}

model_name = 'paci'

if model_name == 'paci':
    scale = 1/1000
else:
    scale = 1

for i, st in enumerate(proto.steps):
    new_seg = []

    if isinstance(st, protocols.VoltageClampStep):
        new_seg.append(f'{st.voltage*scale}')
    else:
        slope = (st.voltage_end - st.voltage_start) / st.duration
        y_int = scale * (st.voltage_start - slope * start_time)

        new_seg.append(f'{slope} * engine.time + {y_int}')

    end_time = (start_time + st.duration * scale)
    new_seg.append(start_time)
    new_seg.append(end_time)

    start_time += st.duration * scale

    all_segments[f'v{i}'] = new_seg

if model_name == 'paci':
    pickle.dump(all_segments, open('./data/opt_proto_s.pkl', 'wb'))
else:
    pickle.dump(all_segments, open('./data/opt_proto_ms.pkl', 'wb'))
