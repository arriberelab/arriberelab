"""
Marcus Viscardi April 1st, 2020

Messing around with multiprocessing for use in assignReadsToGenesX.py script(s)

Major test will be to compare the speeds of writing a bunch of csv lines to a file with:
    1. Multiprocessing, with each process writing to the csv on its own
    2. Multiprocessing, with each process writing to a list/dict/df, then that being writen to csv all at once
    3. Vectorizing (word? not sure), having the process act on the entire df to try and take advantage of the
        faster C backend of vectorized methods

TODO: April 8th, 2020:
    Ideas 1 and 2 will face similar problems of having to synchronize each method step, this could be done using
    something like multiprocessing.Lock:
    (https://docs.python.org/2/library/multiprocessing.html#synchronization-between-processes)
    But I am currently worried that this will overcomplicate the system. I am first going to try to chase after a
    combination of idea 3 and a multiprocessing step.
    The general idea will be to split the df into the number of cores we can utilize. Have each core act on one of
    these groups, then bring the groups back together to be written to the file using a C backend method such
    as pandas.Dataframe.to_csv:
    (https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html)

Main plan will be to implement a toy method for each of the above then run it at a very high repetition
"""
import csv
import timeit
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def timed_csv_write(file, iterations, group_write=True):
    """Comparing CSV writing per line or as one large group"""
    iter_list = range(iterations)
    iter_dict = {str(n): n for n in iter_list}

    def individual_writes():
        with open(file, 'w') as f:
            fwriter = csv.writer(f, delimiter='\t')
            for _ in iter_dict:
                fwriter.writerow(_)

    def group_writes():
        with open(file, 'w') as f:
            fwriter = csv.writer(f, delimiter='\t')
            fwriter.writerows(iter_dict)

    start = timeit.default_timer()
    if not group_write:
        individual_writes()
    else:
        group_writes()
    end = timeit.default_timer()
    if not group_write:
        print(f'INDIVIDUAL:\n\t{end - start} seconds for {iterations} iterations')
    else:
        print(f'GROUP:\n\t{end - start} seconds for {iterations} iterations')

    return end - start


def timed_csv_write_multirun_and_plot():
    """Use to timed_csv_write function from above to run and plot a few different table sizes"""
    iterations = range(7)
    lines_to_write = [10 ** i for i in iterations]
    group_write_time = []
    individual_write_time = []
    file_to_write = 'blah.csv'

    for i in iterations:
        group_write_time.append(timed_csv_write(file_to_write, 10 ** i, group_write=True))
        individual_write_time.append(timed_csv_write(file_to_write, 10 ** i, group_write=False))

    group_per_line = [group_write_time[i] / lines_to_write[i] for i in range(len(iterations))]
    print("Group write:", group_per_line, '\n')
    individual_per_line = [individual_write_time[i] / lines_to_write[i] for i in range(len(iterations))]
    print("Individual Write:", individual_per_line, '\n')

    group_by_ind_ratios = [group_per_line[i] / individual_per_line[i] for i in range(len(iterations))]
    print("Ratios (group/individual):", group_by_ind_ratios)

    plt.plot(iterations, group_write_time, individual_write_time)
    plt.legend(labels=["Group Writing", "Per Line Writing"])
    plt.yscale("linear")
    plt.xscale("linear")
    plt.xlabel("10^i Lines written")
    plt.ylabel("Seconds to write")
    plt.show()


# This will be the major tool to parallelize (in a simple way)
def parallelize_dataframe(df, func, n_cores=4):
    """Parallelize a function call on a dataframe"""
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def vectorized_processing():
    pass


"""
def multiprocessing_with_csv():
    def single_process_to_csv():
        pass
    pass


def multiprocessing_with_list():
    def single_process_to_list():
        pass
    pass
"""


def main(number_of_repititions: int):
    pass


if __name__ == '__main__':
    timed_csv_write('blah.csv', 10**3, group_write=False)
    # timed_csv_write_multirun_and_plot()
