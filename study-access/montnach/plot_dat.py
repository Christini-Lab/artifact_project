import matplotlib.pyplot as plt
import numpy as np


fig, axs = plt.subplots(4, 1, sharex=True)
files = ['INabad_0.600_Rs5.0_Cm20.0']

for fi in files:
    with open(fi) as f:
        lines = f.readlines()

    rows_as_nums = []
    for row in lines:
        rows_as_nums.append([float(n.strip()) for n in row.split('\t') if n != '\n'])
        
    rows_as_nums = np.array(rows_as_nums)


    num_cols = rows_as_nums.shape[1]


    for i in range(0, num_cols-2):
        axs[i].plot(rows_as_nums[:, 0], rows_as_nums[:, i+1])

    print(fi)
        

axs[0].set_ylabel('Vcmd')
axs[1].set_ylabel('Vm')
axs[2].set_ylabel('Im')
axs[3].set_ylabel('Is')

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

axs[3].set_xlabel('Time (ms)')

plt.show()
