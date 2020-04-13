"""
parseFastqToDataframe.py
Marcus Viscardi     April 11, 2020

Trying to parse a fastq file into a pandas dataframe, with the goal of being able to utilize read data faster

The initial goal will be to parse a fastq file into a dataframe with the read info line as one column, the sequence as
    another column, and the confidence as a third column
Later further parsing could allow for more columns from the read info line
Further down the line it would be even better to be able to parse SAM files into a dataframe, which would (maybe) allow
    for acceleration of the assignReadsToGenes functionality

April 11, 2020:
    Success with utilizing pandas dataframe to parse and hold a fastq file. Had no issue paring a fastq of
    14,533,400 lines with my personal desktop, utilizing 1.3GB of RAM once stored as a dataframe. When attempting to
    read in 145,334,000 lines, the memory (16GB) maxed out - but this didn't crash the script. Pandas.read_csv() has a
    memory_low functionality that seems to be acting to somehow avoid crashing. Watching TaskManager, I could see RAM
    utilization bouncing around between 15GB and 15.9GB, out of 16GB total.

"""

import pandas as pd
import argparse
import os, sys

# Pandas default would cut off basically all the columns so:
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 300)


def parseArgs():
    parser = argparse.ArgumentParser(description="Parse a fastq file into a pandas dataframe")
    parser.add_argument('filename', metavar='filename', type=str, help="Path to .fastq file")
    
    args = parser.parse_args()
    file = args.filename
    
    if not os.path.isfile(file):
        print("This is not a file path\n\n\nTERMINATING\n\n")
        sys.exit()
    else:
        print(f"Parsing file @: {file}")
        
    print("Given Arguments:")
    for arg in vars(args):
        print('\t', arg, '=', vars(args)[arg])
    arg_dict = {arg: vars(args)[arg] for arg in vars(args)}
    print(arg_dict)
    return file


def parseFastqToDataframe(file):
    stacked_series = pd.read_csv(file,
                                 sep="\n",
                                 header=None,
                                 squeeze=True,
                                 names=["Values"])
    read_info = stacked_series.iloc[0::4]
    sequence = stacked_series.iloc[1::4]
    confidence = stacked_series.iloc[3::4]
    
    df = pd.concat([read_info.reset_index().Values,
                    sequence.reset_index().Values,
                    confidence.reset_index().Values],
                   axis=1, ignore_index=True)
    df = df.rename(columns={0: "Read Info", 1: "Sequence", 2: "Confidence"})
    print(df)
    print(df.info(memory_usage='deep'))
    return df


if __name__ == '__main__':
    try:
        file_path = parseArgs()
    # If no file is given it will run on the example yeastTestRun in the github directory
    except:
        file_path = "../../../yeastTestRun/SRR4050180TestFile.fastq"
    df = parseFastqToDataframe(file_path)
