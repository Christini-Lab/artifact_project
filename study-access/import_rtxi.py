import h5py
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt


def get_exp_as_df(h5_file, trial_number, cm, t_range=None, trial_type=None):
    """I was going to save the time, voltage and current as a csv,
    but decided not to, because there can be >3million points in 
    the h5 dataset. If you want to make comparisons between trials or
    experiments, call this multiple times.
    """
    data_h5 = h5py.File(h5_file, 'r')

    cm *= 1E-12
    current, voltage = get_current_and_voltage(data_h5, trial_number,
            trial_type=trial_type)

    t_data = get_time_data(data_h5, trial_number)

    d_as_frame = pd.DataFrame({'Time (s)': t_data,
                               'Voltage (V)': voltage,
                               'Current (pA/pF)': current / cm})

    if t_range is not None:
        idx_start = (d_as_frame['Time (s)']-t_range[0]).abs().idxmin()
        idx_end = (d_as_frame['Time (s)']-t_range[1]).abs().idxmin()
        d_as_frame = d_as_frame.copy().iloc[idx_start:idx_end, :]
    

    return d_as_frame


def get_current_and_voltage(f, trial, trial_type=None):

    channels = f[f'Trial{trial}']['Synchronous Data'].keys()

    v_channel = None
    i_channel = None

    for channel in channels:
        if trial_type is not None:
            if trial_type == 'Current Clamp':
                if (('Current Output A' in channel)):
                    i_channel = int(channel.split()[0]) - 1
                if (('Voltage Input V' in channel)):
                    v_channel = int(channel.split()[0]) - 1
                #if (('Analog Output' in channel)):
                #    v_channel = int(channel.split()[0]) - 1
                #if (('Analog Input' in channel)):
                #    i_channel = int(channel.split()[0]) - 1
        else:
            if (('Analog Output' in channel)):
                v_channel = int(channel.split()[0]) - 1
            if (('Analog Input' in channel)):
                i_channel = int(channel.split()[0]) - 1

    if v_channel is None:
        if trial_type is not None:
            for channel in channels:
                if trial_type == 'Current Clamp':
                    if (('Analog Output' in channel)):
                        i_channel = int(channel.split()[0]) - 1
                    if (('Analog Input' in channel)):
                        v_channel = int(channel.split()[0]) - 1

    ch_data = f[f'Trial{trial}']['Synchronous Data']['Channel Data'][()]

    if trial_type is not None:
        if trial_type == 'Current Clamp':
            voltage = ch_data[:, v_channel]
            current = -ch_data[:, i_channel]
            return current, voltage

    channel_1 = ch_data[:, v_channel]
    channel_2 = ch_data[:, i_channel]

    channel_1_test = channel_1[np.logical_not(np.isnan(channel_1))]
    channel_2_test = channel_2[np.logical_not(np.isnan(channel_2))]

    if np.abs(channel_1_test).mean() == 0:
        current = channel_1
        voltage = channel_2

    if np.abs(channel_1_test).std() < np.abs(channel_2_test).mean():
        current = channel_1
        voltage = channel_2
    else:
        current = channel_2
        voltage = channel_1

    avg_early_voltage = voltage[10:100].mean()
    is_voltage_clamp = False 

    if (avg_early_voltage < -.079) and (avg_early_voltage > -.081):
        is_voltage_clamp = True

    if (avg_early_voltage == 0):
        #For funny current
        is_voltage_clamp = True

    if not is_voltage_clamp:
        current = -current

    return current, voltage


def get_time_data(data_h5, trial_number):
    total_time, period = get_time_and_period(data_h5, trial_number)
    ch_data = extract_channel_data(data_h5, trial_number)

    time_array = np.arange(0, len(ch_data[:,0])) * period

    return time_array


def get_time_and_period(data_h5, trial_number):
    start_time, end_time = start_end_time(data_h5, trial_number)
    trial_str = f'Trial{trial_number}'
    total_time = (end_time - start_time) / 1E9
    period = data_h5[trial_str]['Period (ns)'][()] / 1E9

    return total_time, period


def extract_channel_data(data_h5, trial_number):
    trial_str = f'Trial{trial_number}'
    data = data_h5[trial_str]['Synchronous Data']['Channel Data'][()]
    return data


def start_end_time(data_h5, trial_number):
    trial_str = f'Trial{trial_number}'
    start_time = data_h5[trial_str]['Timestamp Start (ns)'][()]
    end_time = data_h5[trial_str]['Timestamp Stop (ns)'][()]
    return start_time, end_time


def main():
    h5_file = './data/h5_files/4_021821_1_alex_control.h5'
    trial_number = 22
    cm = 26

    h5_dat = get_exp_as_df(h5_file, trial_number, cm)


if __name__ == '__main__':
    main()
