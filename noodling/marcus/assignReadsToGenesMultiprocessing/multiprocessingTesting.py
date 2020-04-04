"""
Marcus Viscardi April 1st, 2020

Messing around with multiprocessing for use in assignReadsToGenesX.py script(s)

Major test will be to compare the speeds of writing a bunch of csv lines to a file with:
    1. Multiprocessing, with each process writing to the csv on its own
    2. Multiprocessing, with each process writing to a list/dict/df, then that being writen to csv all at once
    3. Vectorizing (word? not sure), having the process act on the entire df to try and take advantage of the
        faster C backend of vectorized methods

Main plan will be to implement a toy method for each of the above then run it at a very high repetition
"""
import timeit, csv
from itertools import repeat
import pandas as pd
import numpy as np
"""
# Imports from assignReadsToGenes4.py >>>
# import sys, os, collections, csv, common, re, time, pickle, copy, linecache
# import pandas
#from logJosh import Tee
csv.field_size_limit(sys.maxsize)
# <<< Imports from assignReadsToGenes4.py
"""


def timed_csv_write(file, iterations, individual=False):
    iter_list = range(iterations)
    iter_dict = {str(n): n for n in iter_list}

    def individual_writes():
        with open(file, 'a') as f:
            fwriter = csv.writer(f, delimiter='\t')
            for _ in iter_dict:
                fwriter.writerow(_)

    def group_writes():
        with open(file, 'a') as f:
            fwriter = csv.writer(f, delimiter='\t')
            fwriter.writerow((f, iter_dict))

    start = timeit.default_timer()
    if individual:
        individual_writes()
    else:
        group_writes()
    end = timeit.default_timer()
    if individual:
        print(f'INDIVIDUAL:\n\t{end - start} seconds for {iterations} iterations')
    else:
        print(f'GROUP:\n\t{end - start} seconds for {iterations} iterations')



def multiprocessing_with_csv():
    def single_process_to_csv():
        pass
    pass


def multiprocessing_with_list():
    def single_process_to_list():
        pass
    pass


def vectorized_processing():
    pass


def main(number_of_repititions: int):
    pass


if __name__ == '__main__':
    iterations = 10**3
    timed_csv_write('blah.csv', iterations, False)
    timed_csv_write('blah.csv', iterations, True)
