import matplotlib.pyplot as plt
import json
from multiprocessing import cpu_count

if __name__ == '__main__':
    with open('./sp_speeds.json', 'r') as sp_file:
        sp_speed_dict = json.load(sp_file)
    with open('./mp_speeds.json', 'r') as mp_file:
        mp_speed_dict = json.load(mp_file)
    mp_x = [int(x) for x in mp_speed_dict.keys()]
    mp_y = [float(y) for y in mp_speed_dict.values()]
    sp_x = [int(x) for x in sp_speed_dict.keys()]
    sp_y = [float(y) for y in sp_speed_dict.values()]
    
    fig, (axs1, axs2) = plt.subplots(1, 2, sharey=True)
    axs1.set(title='Multiprocessed:\nOne process per DF',
             ylabel=f'Time to process $10^{6}$ DF lines (seconds)',
             xlabel=f'Number of DFs used to subdivide $10^{6}$ lines')
    axs2.set(title='Iterative:\nOne loop per DF',
             xlabel=f'Number of DFs used to subdivide $10^{6}$ lines')
    axs1.axvline(cpu_count(), 0, 1, c='gray', ls='--', label=f'{cpu_count()} '
                                                             f'core processor\n(ideal number of processes?)')
    axs1.legend()
    axs1.plot(mp_x, mp_y, '.-')
    axs2.plot(sp_x, sp_y, '.-')
    plt.show()
