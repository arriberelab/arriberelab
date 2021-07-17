"""
Joshua Arribere, Jan 4, 2018
Converted to python3: April 15, 2020
Marcus Viscardi, July 16, 2021
    Added plotly support, so that you can just highlight genes to see differences

Script to plot select columns of a geneCt file, as output
    from median_normalizer

Input: inFile.geneCt - tab-delimited file of format:
    Gene\tCol1\tCol2...
    gene\tct\tct\t...
    colNames - space-delimited list of column names you
        want to plot. By default, the first column name
        will be the y-axis for all plots

Output: scatter plot with myFavoriteGenes highlighted

run as python3 plotGeneCts_plotly.py inFile.geneCt outPrefix colName1 colName2
"""
import sys, common, csv
from logJosh import Tee


def pdParseGeneCtFile(inFile):
    import pandas as pd
    dataframe = pd.read_csv(inFile, sep="\t")
    dataframe.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
    gene_name_df = pd.read_csv("../geneNames_and_WBGenes.tsv", sep="\t")
    dataframe = dataframe.merge(gene_name_df[["gene_id", "gene_name"]], how="left")
    dataframe["identity"] = dataframe["gene_id"] + " (" + dataframe["gene_name"] + ")"
    return dataframe


def plotlyMkScatterPlot(geneCtDF, columns_to_plot):
    import plotly.express as px
    fig = px.scatter(geneCtDF,
                     x=columns_to_plot[0], y=columns_to_plot[1],
                     hover_name="identity")
    # fig.update_traces(marker=dict(size=12, opacity=0.5, color='#696969',
    #                               line=dict(width=1,
    #                                         color='DarkSlateGrey')
    #                               ),
    #                   selector=dict(mode='markers'))
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    fig.show()


def plotlyMkScatterPlots(geneCtDF, columns_to_plot):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    
    if len(columns_to_plot) == 2:
        plotlyMkScatterPlot(geneCtDF, columns_to_plot)
    else:
        fig = make_subplots(rows=1, cols=len(columns_to_plot)-1)
        
        y_key = columns_to_plot[0]
        for plot_num, x_key in enumerate(columns_to_plot[1:]):
            plot_num += 1
            fig.add_trace(go.Scatter(x=geneCtDF[x_key], y=geneCtDF[y_key],
                                     name=f"{x_key} vs {y_key}",
                                     mode="markers",
                                     hovertext=geneCtDF["identity"],
                                     hoverinfo="text+x+y"),
                          row=1, col=plot_num)
            fig.update_xaxes(title_text=x_key,
                             type="log",
                             row=1, col=plot_num)
            fig.update_yaxes(title_text=y_key,
                             type="log",
                             row=1, col=plot_num)
        fig.update_layout(title_text="Customizing Subplot Axes",
                          showlegend=False)
        fig.show()


def main_plotly_and_pandas(args):
    inFile,outPrefix=args[0:2]
    columns_to_plot = args[2:]
    df = pdParseGeneCtFile(inFile)
    print(df.info())
    plotlyMkScatterPlots(df, columns_to_plot)


if __name__=='__main__':
    Tee()
    main_plotly_and_pandas(sys.argv[1:])
