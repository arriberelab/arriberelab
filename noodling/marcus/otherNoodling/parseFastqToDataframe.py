"""
parseFastqToDataframe.py
Marcus Viscardi     April 11, 2020

Trying to parse a fastq file into a pandas dataframe, with the goal of being able to utilize read data faster

The initial goal will be to parse a fastq file into a dataframe with the read info line as one column, the sequence as
    another column, and the confidence as a third column
Later further parsing could allow for more columns from the read info line
Further down the line it would be even better to be able to parse SAM files into a dataframe, which would (maybe) allow
    for acceleration of the assignReadsToGenes functionality
"""

import pandas as pd
import argparse
import os, sys

parser = argparse.ArgumentParser(description="Parse a fastq file into a pandas dataframe")
parser.add_argument('filename', metavar='filename', type=str, help="Path to .fastq file")


if __name__ == '__main__':
    args = parser.parse_args()
    file = args.filename

    if not os.path.isfile(file):
        print("This is not a file path")
        sys.exit()
    else:
        print(f"Parsing file @: {file}")
