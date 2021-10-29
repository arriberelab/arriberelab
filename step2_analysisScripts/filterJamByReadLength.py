"""
filterJamByReadLength.py
Marcus Viscardi,    October 29, 2021

A script to take a jam and return a jam w/ only reads of a given
length or length range

Run as:
python3 input.jam output.jam read_len1 [read_len2]
    read_len1: required, if read_len2 not passed only reads of this length will survive
    read_len2: optional, if passed then reads from read_len1 to read_len2 will survive
                         including reads of len1 and len2
"""
import sys
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

def main(args):
    print(args)
    in_file, out_file = args[0:2]
    read_len1 = int(args[2])
    jam_df = pd.read_csv(in_file, sep="\t")
    print(f"Original JAM:      {jam_df.shape[0]:>8} reads")
    if len(args) == 4:
        read_len2 = int(args[3])
        filtered_df = jam_df[read_len1 <= jam_df.map_read_seq.apply(len)]
        filtered_df = filtered_df[filtered_df.map_read_seq.apply(len) <= read_len2]
        print(f"Between {read_len1:>2}nt-{read_len2:>2}nt: {filtered_df.shape[0]:>8} reads")
    else:
        filtered_df = jam_df[jam_df.map_read_seq.apply(len) == read_len1]
        print(f"Matching {read_len1:>2}nt:    {filtered_df.shape[0]:>8} reads")
    filtered_df.to_csv(out_file, sep="\t", index=False)
    print("Done!")


if __name__ == '__main__':
    main(sys.argv[1:])
