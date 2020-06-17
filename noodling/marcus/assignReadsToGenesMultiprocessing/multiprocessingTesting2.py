"""
multiprocessingTesting2.py
Marcus Viscardi     April 18, 2020

Revisiting multiprocessing options after 'finishing' with assignReadsToGenesDF.py. The goal will be to test out a way to
    utilize multiprocessing and particularly multiprocessing.Manager().dict() to provide a shared dictionary which will
    contain the chr_key:chr_dataframe values. This will be used to pass one chromosome to each subprocess. The larger
    issue will lie in the writing to this dictionary after functions have been applied. Currently assignReadsToGenes.py
    uses a lot of writes up to the parent dictionary while working in a "for chr_key, chr_df in chr_df_dict.items():"
    loop. This will definitely cause issues if not rewritten for multiprocessing.

First goal will be to get a toy version of the chr_df_dict working with some sort of multiprocessing. A likely option
    will be to take in dataframes with the original dictionary. And write values to a separate dictionary. This will
    provide the advantage of limiting issues with processes reading and writing at the same time, but may cause a
    greater issue of memory over-usage.
"""

import sys
import os
import multiprocessing as mp
import timeit

import matplotlib.pyplot as plt
import json
from multiprocessing import cpu_count

import numpy as np
from pandas import set_option, read_csv, DataFrame, concat, Series
# Pandas default would cut off any columns beyond 5 so:
set_option('display.max_rows', 50)
set_option('display.max_columns', 20)
set_option('display.width', 300)


def buildFake_ChrDfDict(chrs: int = 10, length: int = 10**5):
    chr_df_dict = {}
    for chr in range(chrs):
        df = DataFrame()
        df.insert(0, 'Numbers', [x[0] for x in np.random.randint(0, 1000000, size=(length//chrs, 1))])
        chr_df_dict[f'Chr-{chr}'] = df
    # print(f"Finished building fake dictionary of chr:chr_df with keys:\n{chr_df_dict.keys()}")
    return chr_df_dict


def df_edit(df, extra_loops_per_line: int = 10) -> DataFrame:
    
    def single_line_edit(value_to_edit: int) -> int:
        for _ in range(extra_loops_per_line):
            _ = value_to_edit // 15
        return value_to_edit // 15
    
    new_df = DataFrame()
    new_df['Numbers'] = Series(df.apply(lambda x: single_line_edit(x['Numbers']), axis=1).tolist(), index=df.index)
    return new_df


def one_by_one(chromosomes=10, length=10**6):
    start = timeit.default_timer()
    chr_df_dict = buildFake_ChrDfDict(chrs=chromosomes, length=length)
    for chr_key, df in chr_df_dict.items():
        # print(f"{chr_key}, example value change:\n\t{chr_df_dict[chr_key].iloc[0]['Numbers']:>6}", end='')
        updated_df = df_edit(df)
        chr_df_dict[chr_key] = updated_df
        # print(f" -> {chr_df_dict[chr_key].iloc[0]['Numbers']:<6}")
    end = timeit.default_timer()
    # print(f'one_by_one: {end - start} seconds')
    return end - start


def workerProcess(chr, df_dict_in, df_dict_out):
    # print(f"Process {os.getpid():>5} working on {chr}")
    new_df = df_edit(df_dict_in[chr])
    df_dict_out[chr] = new_df


def all_at_once(chromosomes=10, length=10**6):
    chr_df_dict = buildFake_ChrDfDict(chrs=chromosomes, length=length)
    start = timeit.default_timer()
    manager = mp.Manager()
    shared_input_dict, shared_output_dict = manager.dict(chr_df_dict), manager.dict()
    process_list = []
    for chr_key, chr_df in shared_input_dict.items():
        process = mp.Process(target=workerProcess, args=[chr_key, shared_input_dict, shared_output_dict])
        process_list.append(process)
    # print(f"Processors assigned, {len(process_list)} at the ready\n"
    #       f"\tOutput dictionary has {len(shared_output_dict.keys())} items")
    for pro in process_list:
        pro.start()
    # print(f"Processes started using process.start()\n\tOutput dictionary has {len(shared_output_dict.keys())} items")
    for pro in process_list:
        # print(f'Joining process {pro}')
        pro.join()
    # print(f"Processes joined using process.join()\n\tOutput dictionary has {len(shared_output_dict.keys())} items")
    end = timeit.default_timer()
    # print(f'all_at_once: {end - start} seconds')
    return end -start


if __name__ == '__main__':
    # one_by_one()
    mp_speed_dict = {}
    sp_speed_dict = {}
    total_lines = 10**6
    max_divisions = 24
    min_divisions = 5
    sp_runs = 2

    # Run {sp_runs} single processor runs to be used for an average
    print(f"Running {sp_runs} single processor run(s) to average as a baseline:")
    for i in range(sp_runs):
        speed = one_by_one(chromosomes=i, length=total_lines)
        sp_speed_dict[i] = speed
        print(i+1, end=' ')
    print()

    print(f"Running multiprocessed runs, using {min_divisions} to {max_divisions} processors:")
    for n in range(max_divisions)[min_divisions:]:
        processes = n+1
        speed = all_at_once(chromosomes=processes, length=total_lines)
        mp_speed_dict[processes] = speed
        print(processes, end=' ')
    print()
    for proc, speed in mp_speed_dict.items():
        print(f'When {total_lines} lines are split into {proc:>2} dataframes of {total_lines//proc} '
              f'lines each (w/ a processor per DF), it took {speed:2.2f} seconds to process them all')
    mp_x = [int(x) for x in mp_speed_dict.keys()]
    mp_y = [float(y) for y in mp_speed_dict.values()]
    sp_x = [int(x) for x in sp_speed_dict.keys()]
    sp_y = [float(y) for y in sp_speed_dict.values()]
    
    # with open('./mp_speeds.json', 'w') as mp_file:
    #     json.dump(mp_speed_dict, mp_file)
    # with open('./sp_speeds.json', 'w') as sp_file:
    #     json.dump(sp_speed_dict, sp_file)
    #
    # with open('./sp_speeds.json', 'r') as sp_file:
    #     sp_speed_dict = json.load(sp_file)
    # with open('./mp_speeds.json', 'r') as mp_file:
    #     mp_speed_dict = json.load(mp_file)

    ####################################################################################################################
    # mp_x = [int(x) for x in mp_speed_dict.keys()]
    # mp_y = [float(y) for y in mp_speed_dict.values()]
    # sp_x = [int(x) for x in sp_speed_dict.keys()]
    # sp_y = [float(y) for y in sp_speed_dict.values()]
    sp_average = sum(sp_y)/len(sp_y)

    fig, axs = plt.subplots(1, 1, sharey=True)
    axs.set(title='Multiprocessed:\nOne process per DF',
             ylabel=f'Time to process $10^{6}$ DF lines (seconds)',
             xlabel=f'Number of DFs used to subdivide $10^{6}$ lines')
    axs.axvline(cpu_count(), 0, 1, c='gray', ls='--', label=f'{cpu_count()} '
                                                            f'core processor\n(ideal number of processes?)')
    axs.axhline(sp_average, 0, 1, c='gray',)
    axs.legend()
    axs.plot(mp_x, mp_y, '.-')
    plt.show()
