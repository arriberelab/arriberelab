#!/usr/bin/env python3

class CommandLine():

    """
    Handle the command line, usage and help requests
    CommandLine uses argparse, now standard in 2.7 and beyond
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die
    attributes:
    all arguments received from the commandline using .add_argument will be
    available within the .args attribute of object instantiated from CommandLine
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    """
    def __init__(self, inOpts=None):
        """
        Implement a parser to interpret the command line argv string using argparse.
            --> Argparse adjusted so that default is minGene = 100
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Finds the largest ORF in a DNA sequence',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store_true', help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), action='store',
                                 help='minimum Gene length', default = 100)
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

def main(inCL=None):
    if inCL is None:
        myCommandLine = CommandLine()
    print(myCommandLine.args)

if __name__ == "__main__":
    main()

