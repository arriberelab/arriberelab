#!/usr/bin/env python3
"""
Aug 27, 2013 Joshua Arribere
Converted to python 3: Mar 24, 2020

Script to perform automatic reporting of another python script when invoked. This script will:
(1) Create a file with filename including script used, and date/time
(2) Will create two versions of the file above, one .log, the other .cmd.
    .log contains stdout
    .cmd contains the command run by the

The goal in doing this is to know the precise origins of any given file. I will at least know the command
that was run to generate a file at the time that file was generated.

This was originally adapted from a script from Jason Merkin.
EDIT: Oct 4, 2013 - JOSH revised to make output logs in the directory where python is called. This makes
    trying to identify the command used to make a given file that much easier.
"""
import sys
import os
import datetime


# from http://mail.python.org/pipermail/python-list/2007-May/438106.html
# and http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python

class Tee(object):
    def __init__(self, filename=None, script_name=os.path.basename(sys.argv[0]),
                 time=None, mode='w'):
        time = time or str(datetime.datetime.now()).replace(' ', '_')
        try:
            os.mkdir('logs')
        except OSError:  # then the log directory already exists
            pass
        filename = filename or os.getcwd() + f'/logs/{script_name}_{time}'
        # filename = filename or f'/home/josh/data/logs/{script_name}_{time}'
        print(f'logging to {filename}')
        self.file = open(f'{filename}.log', mode)
        self.stdout = sys.stdout
        sys.stdout = self
        with open(f'{filename}.cmd', 'w') as f:
            f.write(f'From directory: {os.getcwd()}\n')  # let's you know where you are
            f.write(' '.join(sys.argv[0:]) + '\n')  # writes the command that was run.

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        pass


if __name__ == '__main__':
    # Tee()
    Tee(script_name=os.path.basename(__file__))

    print('testing')
    print(1, 2)
    print('passed')
