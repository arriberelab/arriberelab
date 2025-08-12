"""
plotlyGeneCts_colorByContinuous.py
Joshua Arribere, Jan 4, 2018
Converted to python3: April 15, 2020
Marcus Viscardi, July 16, 2021
    Added plotly support, so that you can hover over genes to see what they are
Marcus Viscardi, Oct 20, 2022
    Copied most of this script over from plotlyGeneCts.py, but adding the ability
    to merge with an input dataframe file and coloring genes based on that!

Script to plot select columns of a geneCt file, as output
    from median_normalizer

Input: inFile.geneCt - tab-delimited file of format:
    Gene\tCol1\tCol2...
    gene\tct\tct\t...
    colNames - space-delimited list of column names you
        want to plot. By default, the first column name
        will be the y-axis for all plots

Output: scatter plot with myFavoriteGenes highlighted, also it will spit out an HTML file that is interactive!

run as (spaces separating arguments):
python3 plotlyGeneCts.py inFile.geneCt
                         outputFile(.pdf/.svg/.png)
                         dataframeFileWContinuousData(.csv/.tsv/.parquet)
                         columnNameWithContinuousData
                         X-axis
                         Y-axis
"""
import sys, os

# Hardcoded path to the tsv file that has gene names & WBGene IDs:
#   This file should get copied over by anyone that pulled this script from github!
PATH_TO_GENE_CONVERTER = f"{os.path.abspath(os.path.dirname(__file__))}/../geneNames_and_WBGenes.tsv"
ADD_GENE_NAME_FLAG = True


def pdParseGeneCtFile(inFile, dataframeFileWithContinuousData, continuousDataColumn):
    import pandas as pd
    dataframe = pd.read_csv(inFile, sep="\t")
    dataframe.rename(columns={dataframe.columns[0]: "gene_id"}, inplace=True)
    
    if dataframeFileWithContinuousData.endswith(".tsv"):
        cont_df = pd.read_table(dataframeFileWithContinuousData)
    elif dataframeFileWithContinuousData.endswith(".csv"):
        cont_df = pd.read_csv(dataframeFileWithContinuousData)
    elif dataframeFileWithContinuousData.endswith(".parquet"):
        cont_df = pd.read_parquet(dataframeFileWithContinuousData)
    else:
        raise NotImplemented(f"Please provide a dataframe file that ends in .tsv/.csv/.parquet,"
                             f"\n\tThe file you provides is: {dataframeFileWithContinuousData}")
    try:
        cont_df = cont_df[['gene_id', continuousDataColumn]]
    except KeyError:
        raise KeyError(f"'gene_id' or '{continuousDataColumn}' not found in the dataframe file!"
                       f"\nThese are required to merge with the geneCts file and plot!!")
    dataframe = dataframe.merge(cont_df, how="left")
    if ADD_GENE_NAME_FLAG:
        try:
            gene_name_df = pd.read_csv(PATH_TO_GENE_CONVERTER, sep="\t")
        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find the gene_id to gene_name converter file at: "
                                    f"{PATH_TO_GENE_CONVERTER}\nThis should have been loaded with the github repo! "
                                    f"But you may need to change the hardcoded path in this script.\n\t"
                                    f"For example, on geneGenie you can change line 31 of this script to say:"
                                    f"\n\t\tPATH_TO_GENE_CONVERTER = "
                                    f"\"/data14/roton/scripts/arriberelab/geneNames_and_WBGenes.tsv\""
                                    f"\n\tAnd this will solve the FileNotFoundError!")
        dataframe = dataframe.merge(gene_name_df[["gene_id", "gene_name"]], how="left")
        dataframe["identity"] = dataframe["gene_name"] + " (" + dataframe["gene_id"] + ")"
    else:
        dataframe["identity"] = dataframe["gene_id"]
    return dataframe


def plotlyMkScatterPlot(geneCtDF, output_file, columns_to_plot, continuousDataColumn):
    import plotly.express as px
    from plotly.colors import DEFAULT_PLOTLY_COLORS as DEFAULT_COLORS
    fig = px.scatter(geneCtDF,
                     x=columns_to_plot[1], y=columns_to_plot[0],
                     color=continuousDataColumn,
                     color_continuous_scale="rainbow",
                     opacity=0.8,
                     # color_discrete_sequence=DEFAULT_COLORS,
                     hover_name="identity")
    fig.update_traces(marker=dict(size=4))
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    fig.write_image(output_file)
    fig.write_html(output_file + ".html")
    fig.show()


def main_plotly_and_pandas(args):
    inFile = args[0]
    outFile = args[1]
    dataframeFileWithContinuousData = args[2]
    continuousDataColumn = args[3]
    columns_to_plot = args[4:]
    df = pdParseGeneCtFile(inFile, dataframeFileWithContinuousData, continuousDataColumn)
    print(df.info())
    plotlyMkScatterPlot(df, outFile, columns_to_plot, continuousDataColumn)


if __name__ == '__main__':
    main_plotly_and_pandas(sys.argv[1:])
