#!/usr/bin/env python3
"""
filterBAMbyJAM.py
Marcus Viscardi,    September 08, 2021

Simple script to load the bam from STAR and a final jam
    file, then output a new filtered bam file.
The reason we need this script is because there are several
    steps of read filtering that we do after the STAR step
    of the pipeline. Because the BAM file is made before
    these steps, it is oblivious to what reads we may have
    dropped. Despite this we use the BAM file for things like
    IGV!

Run:
    python3 filterBAMbyJAM.py path/to/bam path/to/jam path/to/new_bam
"""
import pandas as pd
import subprocess as sub
import io
import sys
from collections import namedtuple

BamHeadersAndDf = namedtuple("BamHeadersAndDf", ["headers", "df"])


def bam_to_df(bam_path) -> BamHeadersAndDf:
    # I did a few checks and found that loading the BAM file
    # like this is the sam as loading the SAM file
    
    o = 'object'
    i = 'int64'
    chr_cat = 'category'
    bam_dtypes_list = [o, i, chr_cat, i, i, o, o, i, i, o, o, o, o, o, o]
    bam_dtypes_dict = {i: bam_dtypes_list[i] for i in range(15)}

    output = sub.check_output(f"samtools view {bam_path}", shell=True)
    df = pd.read_csv(io.BytesIO(output),
                     encoding='utf8',
                     sep="\t",
                     names=range(15),
                     dtype=bam_dtypes_dict,
                     )
    df = df.rename(columns={0: 'read_id',
                            1: 'strand',
                            2: 'chr',
                            3: 'chr_pos',
                            4: 'mapq',
                            5: 'cigar',
                            6: 'rnext',
                            7: 'pnext',
                            8: 'tlen',
                            9: 'read_seq',
                            10: 'phred_qual',
                            11: 'NH',
                            12: 'HI',
                            13: 'AS',
                            14: 'nM'})
    
    header = sub.check_output(f"samtools view -H {bam_path}", shell=True).decode("utf-8")
    output = BamHeadersAndDf(header, df)
    return output


def jam_to_df(jam_path) -> pd.DataFrame:
    df = pd.read_csv(jam_path, sep="\t")
    return df


def filter_bam_with_jam(bam_obj: BamHeadersAndDf, jam_df: pd.DataFrame) -> BamHeadersAndDf:
    filtered_bam_df = bam_obj.df[bam_obj.df.read_id.isin(jam_df.read_id)].reset_index(drop=True)
    retained_headers = bam_obj.headers
    return BamHeadersAndDf(retained_headers, filtered_bam_df)


def save_bam_obj(bam_obj: BamHeadersAndDf, output_path: str):
    header, df = bam_obj
    buffer = header + df.to_csv(sep="\t",
                                header=False,
                                index=False)
    sub.run(f"samtools view -S -b - > {output_path}", input=buffer.encode('utf-8'), shell=True)


def main(bam_path, jam_path, output_path):
    print(f"\nStarted loading BAM file. . .", end="\r")
    original_bam = bam_to_df(bam_path)
    print(f"Successfully loaded original BAM file from: {bam_path}\n"
          f"{original_bam.df.shape[0]:>15} reads found\n")
    print(f"Started loading JAM file. . .", end="\r")
    jam_df = jam_to_df(jam_path)
    print(f"Successfully loaded JAM file from: {jam_path}\n"
          f"{jam_df.shape[0]:>15} reads found\n")
    print(f"Started filtering reads. . .", end="\r")
    new_bam = filter_bam_with_jam(original_bam, jam_df)
    print(f"Finished filtering reads based off content of JAM file\n"
          f"{new_bam.df.shape[0]:>15} reads passed\n")
    print(f"Started saving new BAM file. . .", end="\r")
    save_bam_obj(new_bam, output_path)
    print(f"Saved reads that passed to new BAM file at: {output_path}\n")


def troubleshoot_run():
    test_bam_path = "/data16/old-users/john/working/200716_JAA-035/" \
                    "200716_JK_0043_N8_processed/" \
                    "200716_JK_0043_15-18.finalMapped.Aligned.out.sorted.bam"
    
    test_jam_path = "/data16/old-users/john/working/200716_JAA-035/" \
                    "200716_JK_0043_N8_processed/" \
                    "200716_JK_0043_15-18.allChrs.jam"
    
    test_output_path = "./deleteme2.bam"
    
    main(test_bam_path, test_jam_path, test_output_path)
    test_finished = bam_to_df(test_output_path)
    print(test_finished.headers)
    print(test_finished.df)


if __name__ == '__main__':
    # troubleshoot_run()
    bam, jam, out = sys.argv[1:]
    main(bam, jam, out)
