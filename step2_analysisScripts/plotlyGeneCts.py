"""
plotlyGeneCts.py
Joshua Arribere, Jan 4, 2018
Converted to python3: April 15, 2020
Marcus Viscardi, July 16, 2021
    Added plotly support, so that you can just hover over genes to see what they are

Script to plot select columns of a geneCt file, as output
    from median_normalizer

Input: inFile.geneCt - tab-delimited file of format:
    Gene\tCol1\tCol2...
    gene\tct\tct\t...
    colNames - space-delimited list of column names you
        want to plot. By default, the first column name
        will be the y-axis for all plots

Output: scatter plot with myFavoriteGenes highlighted, also it will spit out an HTML file that is interactive!

run as:
python3 plotlyGeneCts.py inFile.geneCt outputFile(.pdf/.svg/.png) myFavWBGenes.txt [space seperated list of cols]
                                    First col will be the Y-axis for all, following will be X-axes^

If you don't want to highlight any genes:
    just write the word NONE instead of a path to a myFavGenes.txt file!
"""
import sys, os

# Hardcoded path to the tsv file that has gene names & WBGene IDs:
#   This file should get copied over by anyone that pulled this script from github!
PATH_TO_GENE_CONVERTER = f"{os.path.abspath(os.path.dirname(__file__))}/../geneNames_and_WBGenes.tsv"
ADD_GENE_NAME_FLAG = True


def pdParseGeneCtFile(inFile, fav_genes_file_path):
    import pandas as pd
    dataframe = pd.read_csv(inFile, sep="\t")
    dataframe.rename(columns={dataframe.columns[0]: "gene_id"}, inplace=True)
    
    if fav_genes_file_path.lower() != 'none':
        # The following loads the list of fav genes
        fav_gene_list = pd.read_csv(fav_genes_file_path, names=['gene_id']).gene_id.to_list()
        # Then this makes a column w/ True if the gene was in the favorites list
        #   and False if not. We can then color everything by that!
        dataframe['fav_gene'] = dataframe['gene_id'].isin(fav_gene_list)
    else:
        dataframe['fav_gene'] = False

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


def plotlyMkScatterPlot(geneCtDF, output_file, columns_to_plot):
    import plotly.express as px
    from plotly.colors import DEFAULT_PLOTLY_COLORS as DEFAULT_COLORS
    fig = px.scatter(geneCtDF,
                     x=columns_to_plot[1], y=columns_to_plot[0],
                     color='fav_gene',
                     color_discrete_sequence=DEFAULT_COLORS,
                     hover_name="identity")
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    fig.write_image(output_file)
    fig.write_html(output_file + ".html")
    fig.show()


def plotlyMkScatterPlots(geneCtDF, output_file, columns_to_plot):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    from plotly.colors import DEFAULT_PLOTLY_COLORS as DEFAULT_COLORS

    if len(columns_to_plot) == 2:
        plotlyMkScatterPlot(geneCtDF, output_file, columns_to_plot)
    else:
        fig = make_subplots(rows=1, cols=len(columns_to_plot)-1)

        y_key = columns_to_plot[0]
        for plot_num, x_key in enumerate(columns_to_plot[1:]):
            plot_num += 1
            fig.add_trace(go.Scatter(x=geneCtDF[x_key], y=geneCtDF[y_key],
                                     marker_color=geneCtDF['fav_gene'].replace({
                                         True: DEFAULT_COLORS[1],
                                         False: DEFAULT_COLORS[0]}),
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
        fig.update_layout(showlegend=False)
        fig.write_image(output_file)
        fig.write_html(output_file + ".html")
        fig.show()


def main_plotly_and_pandas(args):
    inFile, outFile = args[0:2]
    fav_genes_path = args[2]
    columns_to_plot = args[3:]
    df = pdParseGeneCtFile(inFile, fav_genes_path)
    print(df.info())
    plotlyMkScatterPlots(df, outFile, columns_to_plot)


if __name__ == '__main__':
    main_plotly_and_pandas(sys.argv[1:])
