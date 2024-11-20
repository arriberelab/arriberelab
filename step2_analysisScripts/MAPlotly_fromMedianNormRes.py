"""
MAPlotly_fromMedianNormRes.py
Marcus Viscardi,    November 20, 2024

Quick script to take the output from medianNormalizerWithDESeq2_jams.py and turn it into a MAPlot!

Inputs:
    <path to med norm output>.DESeqgeneCts
    <output prefix for html/svg/png>
"""

import argparse
import plotly.express as px
import pandas as pd
from pathlib import Path

def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("medNormRes", type=Path,
                        help="File containing the median normalized gene counts from DESeq2.")
    parser.add_argument("outPrefix", type=str,
                        help="Prefix for the output plot files.")
    parser.add_argument("-t", "--html", action="store_true",
                        help="If set, will output an html file.")
    parser.add_argument("-s", "--svg", action="store_true",
                        help="If set, will output an svg file.")
    parser.add_argument("-p", "--png", action="store_true",
                        help="If set, will output a png file.")
    parser.add_argument("--no-save", action="store_true",
                        help="If set, will not save any files.")
    
    parser.add_argument("-v", "--p-val", type=float, default=0.05,
                       help="P-value cutoff for significance in the MA plot.")
    parser.add_argument("--use-unadj", action="store_true",
                        help="If set, will use the unadjusted p-value for the MA plot.")
    parser.add_argument("-m", "--min-avg-count", type=float, default=1,
                        help="Minimum average count for a gene to be included in the plot.")
    return parser.parse_args()


def make_MAPlot(medNormRes: Path, outPrefix: str,
                p_value_cutoff: float,
                no_save: bool = False,
                use_unadj: bool = False,
                html: bool = False,
                svg: bool = False,
                png: bool = False,
                min_avg_count: float = 1):
    # Read in the data
    medNormDF = pd.read_csv(medNormRes)
    medNormDF.rename(columns={"Unnamed: 0": "Gene"}, inplace=True)
    if use_unadj:
        sig_col_name = f'significant (p-val > {p_value_cutoff})'
        p_col = 'pvalue'
    else:
        sig_col_name = f'significant (padj > {p_value_cutoff})'
        p_col = 'padj'
    medNormDF[sig_col_name] = medNormDF[p_col] <= p_value_cutoff
    
    medNormDF = medNormDF[medNormDF['baseMean'] >= min_avg_count]
    
    # Make the MA plot
    fig = px.scatter(medNormDF,
                     x="baseMean",
                     y="log2FoldChange",
                     color=sig_col_name,
                     color_discrete_map={True: "red", False: "black"},
                     hover_name="Gene",
                     hover_data=["padj",
                                 "pvalue",],
                     log_x=True,
                     title=f"MAPlot for {medNormRes}")
    # Save the plot
    if not no_save:
        write_all = True
        if html:
            fig.write_html(outPrefix + ".html")
            write_all = False
        if svg:
            fig.write_image(outPrefix + ".svg")
            write_all = False
        if png:
            fig.write_image(outPrefix + ".png")
            write_all = False
        if write_all:
            fig.write_html(outPrefix + ".html")
            fig.write_image(outPrefix + ".svg")
            fig.write_image(outPrefix + ".png")
    fig.show()

def main():
    args = get_args()
    
    make_MAPlot(args.medNormRes, args.outPrefix,
                args.p_val, args.no_save, args.use_unadj,
                args.html, args.svg, args.png,
                args.min_avg_count)
    print("Done!")

if __name__ == '__main__':
    main()
