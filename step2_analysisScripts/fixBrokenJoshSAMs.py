"""
fixBrokenJoshSAMs.py
Marcus Viscardi     April 26, 2021

This script should be unnecessary after the update I pushed earlier today,
    but it will be useful for fixing old files that came out of pipeline9
    
This script edits the gene and gene_string columns of the joshSAM files so
    that they match those of the old versions. Hopefully this will mean that
    older scripts will be easy to use with joshSAMs coming out of the newer
    pipeline.
    
THIS DOES NOT MEAN YOU SHOULD WRITE NEWER SCRIPTS TO USE joshSAMS FILES!!
    The jam files retain a lot more information that could be useful down
    the line. So please do you're best to write anything new to work with
    the jam files. <3
"""
from pandas import read_csv
from csv import QUOTE_NONE
from argparse import ArgumentParser


if __name__ == '__main__':
    parser = ArgumentParser(description='Wrapper for the handling of libraries starting from fastq files.')
    # Required Arguments:
    parser.add_argument('joshSAM_input', metavar='joshSAM_input',
                        type=str, help='joshSAM file from previous version of pipeline9 that were broken')
    parser.add_argument('output_file', metavar='output_file', type=str,
                        help='Name for the output joshSAM file (don\'t include the .joshSAM)')
    args = parser.parse_args()

    # Load your file:
    df = read_csv(args.joshSAM_input, sep="\t")
    print(f"Loaded joshSAM file from: {args.joshSAM_input}")
    # Add the :ASorS and replace pipe separators with tabs:
    df["gene_string"] = df.apply(lambda row: row["gene_string"].replace("|", f":{row['gene'].split(':')[-1]}\t") + f":{row['gene'].split(':')[-1]}", axis=1)
    # Remove the :ASorS from the gene IDs:
    df["gene"] = df.apply(lambda row: row["gene"].split(":")[0], axis=1)
    
    # Save the file!!
    df.to_csv(f"{args.output_file}.joshSAM", sep="\t",
              header=False, index=False, quoting=QUOTE_NONE, quotechar="",  escapechar="\\")
    print(f"Done! Look for your fixed joshSAM file at: {args.output_file}.joshSAM")
