"""
argParseTesting.py
Marcus Viscardi,    October 04, 2021

Trying to see if I can use a single flag to add the 4 variables
    needed for the QC graphic!
"""
import argparse

# Absolute defaults are overwritten by the given settings file and any command line arguments given
ABSOLUTE_DEFAULT_DICT = {'keepNonUnique': False, 'outputJoshSAM': False,
                         'printArgs': False, 'filterMap': False, 'regenerate': False,
                         'cores': 7, 'misMatchMax': 0,
                         'optString': '--outFilterScoreMinOverLread 1 '
                                      '--outFilterMatchNminOverLread 1 '
                                      '--outReadsUnmapped Fastx '
                                      '--outSJfilterOverhangMin 6 6 6 6',
                         'misMatchMax2': 3,
                         'optString2': f'--outFilterScoreMin 14 '
                                       f'--outFilterScoreMinOverLread 0.3 '
                                       f'--outFilterMatchNmin 14 '
                                       f'--outFilterMatchNminOverLread 0.3 '
                                       f'--outReadsUnmapped Fastx '
                                       f'--outSJfilterOverhangMin 1000 1000 1000 1000 ',
                         'genomeDir2': False, 'genomeAnnots2': False,
                         'heatmapWindows': [-21, 21, -30, 12]}


# Stealing below from the main pipelineWrapper:
def parseArguments():
    """
    Outputs a dictionary of argument keys (taken either from the metavar parameter below or the --name for flags)
        and the passed argument's item. This can be used with
    """
    parser = argparse.ArgumentParser(description='Wrapper for the handling of libraries starting from fastq files.')
    # Required Arguments:
    parser.add_argument('settings', metavar='settings', type=str,
                        help='A text file containing a line-delimited input of adaptorSeq, genomeDir, and genomeAnnots.'
                             ' Plus any of the following optional arguments. All in arg|value format')
    # Optional Arguments: (None default here allows us to not pass anything that isn't given by user.
    #                      This helps to simplify settings.txt/default overwrites further down the line.)
    parser.add_argument('--umi5', metavar='umi5', type=int, default=None,
                        help='The number of nucleotides to be trimmed from the 5prime end of reads.')
    parser.add_argument('--umi3', metavar='umi3', type=int, default=None,
                        help='The number of nucleotides to be trimmed from the 3prime end of reads.')
    parser.add_argument('--minimumReadLength', '--min', type=int, default=None,
                        help='The shortest reads to be mapped to the genome.')
    parser.add_argument('--maximumReadLength', '--max', type=int, default=None,
                        help='The longest reads to be mapped to the genome.')
    parser.add_argument('--adaptorSeq', '--adaptor', metavar='adaptor', type=str, default=None,
                        help='3\' adaptor to be trimmed off.')
    parser.add_argument('--genomeDir', metavar='genomeDir', type=str, default=None,
                        help='Genome directory where STAR index can be found.')
    parser.add_argument('--genomeAnnots', metavar='genomeAnnots', type=str, default=None,
                        help='Genome annotations (gtf format).')
    parser.add_argument('--cores', metavar='cores', type=int, default=None,
                        help='Number of cores to use.')
    parser.add_argument('--misMatchMax', metavar='misMatchMax', type=int, default=None,
                        help='Number of mismatches to tolerate during mapping.')
    parser.add_argument('--heatmapWindows', nargs=4, metavar=('upStart', 'downStart', 'upStop', 'downStop'),
                        help="Change the window sizes of the heatmaps that come out of the pipeline. This"
                             "argument takes 4 values seperated by spaces: upStart downStart upStop downStop")
    # Flag Arguments: (just add these as tags to change pipeline functionality)
    parser.add_argument('-k', '--keepNonUnique', action='store_true',
                        help="Boolean flag to allow non-unique reads in the final .jam file(s).")
    parser.add_argument('-j', '--outputJoshSAM', action='store_true',
                        help="Boolean flag to also output joshSAM format files in addition to the jam")
    parser.add_argument('-p', '--printArgs', action='store_true',
                        help="Boolean flag to show how arguments are overwritten/accepted")
    parser.add_argument('-f', '--filterMap', action='store_true',
                        help="Boolean flag to perform filter mapping")
    parser.add_argument('-r', '--regenerate', action='store_true',
                        help="Boolean flag to ignore previously produced files and generate all files anew")

    # Spit out namespace object from argParse
    args = parser.parse_args()
    # Quickly convert Namespace object to dictionary
    arg_dict = {arg: vars(args)[arg] for arg in vars(args)}

    if arg_dict['printArgs']:
        # Print arguments, specifies arguments which will not be passed in the arg_dict
        print("\nGiven Arguments (ArgParse):")
        for key, arg in arg_dict.items():
            if not arg:
                print(f"\t{key} = {arg} -> (Will not be passed)")
            else:
                print(f"\t{key} = {arg}")
    # Recreate dict without arguments that did not receive any input
    arg_dict = {k: v for k, v in arg_dict.items() if v is not None and v}
    return arg_dict


def parseSettings(settings, printArgs=False, **other_kwargs):
    """
    Will loop through and replace variables that are None
    """
    # first, parse the settings file to a dictionary called settingsDict
    settingsDict = {}
    with open(settings, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip()
                if line != '':
                    line = line.split('|')
                    if len(line) == 2:
                        settingsDict[line[0]] = line[1]
                    else:
                        print("\033[31;1m\nRemove pipes ('|') from settings "
                              "file arguments (or rewrite parser)\n\033[0m")
                        exit()
    if printArgs:
        print(f"\nSettings Arguments (file: '{settings}')")
        for key, arg in settingsDict.items():
            if not arg:
                print(f"\t{key} = {arg} -> (Will not be passed)")
            else:
                print(f"\t{key} = {arg}")
    settingsDict = {k: v for k, v in settingsDict.items() if v is not None and v is not ''}
    return settingsDict


def combineSettingsAndArguments():
    absoluteDefDict = ABSOLUTE_DEFAULT_DICT
    argDict = parseArguments()
    settingsDict = parseSettings(**argDict)
    finalArgDict = {}

    finalArgDict.update(absoluteDefDict)
    finalArgDict.update(settingsDict)
    finalArgDict.update(argDict)
    print("\033[1m\nPipeline Arguments:")
    # Absolute defaults overwritten by settings.txt then overwritten by CLI args
    for key, arg in finalArgDict.items():
        print(f"\t{key} = {arg}")
        if arg == "True":
            finalArgDict[key] = True
        elif arg == "False":
            finalArgDict[key] = False
        elif key == "heatmapWindows":
            if isinstance(arg, list):
                finalArgDict[key] = list(map(int, arg))
            else:
                finalArgDict[key] = list(map(int, arg.split()))
        else:
            try:
                finalArgDict[key] = int(arg)
            except ValueError:
                finalArgDict[key] = str(arg)
            except TypeError:
                finalArgDict[key] = str(arg)
    # print(f"cutadapt version:")
    # os.system('cutadapt --version')
    # print(f"STAR version:")
    # os.system('STAR --version')
    print('\n\033[0m\n')
    return finalArgDict


if __name__ == '__main__':
    argument_dict = combineSettingsAndArguments()
    print(argument_dict)
