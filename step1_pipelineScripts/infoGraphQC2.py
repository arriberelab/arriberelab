"""
Joshua Arribere, April 16, 2020

Script to make an infograph summarizing features of a library.

Input: file.jam - josh sam file
    lower, upper - shortest and longest read length to consider, inclusive
    upStart,downStart - will look this many nts up/downsteram of start
        codon
    upStop,downStop - will look this many nts up/downstream of stop codon

Output: a file showing the metaStart, metaStop, fraction of reads of each
    read length in [lower,upper], and the fraction in frame for each read
    length

run as python infoGraphQC2.py file.jam lower upper upStart downStart upStop
    downStop outPrefix
"""
import sys, common
import metaStartStopHeatMap2, readLengthAndPhasingAnalysis2
import seaborn, matplotlib, pandas
from logJosh import Tee

def mkPlot(dfStart,dfStop,df,outPrefix):
    """
    Will make heatmaps of dfStart and dfStop, then bar graphs of df, which
    contains read lengths and frame information.
    """
    ##first initialize the style
    seaborn.set_style("whitegrid", {'grid.linestyle': '--'})
    ##initialize the figures, axes
    fig,axs=matplotlib.pyplot.subplots(nrows=4,ncols=1,
        figsize=(6.4,4.8*2))
    ##set some parameters for the heatmap colorkey
    cbar_ax = fig.add_axes([.86,  # x origin (was .81)
                            .55,  # y origin
                            .03,  # x width
                            .3],  # y width
                           title='RPM')
    ##plot the metaStart
    seaborn.heatmap(dfStart,ax=axs[0],cmap="YlGnBu",
        square=True,
        linewidths=.5,vmax=400,
        xticklabels=3,#makes the labels every 3 nts
        cbar=True,cbar_ax=cbar_ax)
    ##plot the metaStop
    seaborn.heatmap(dfStop,ax=axs[1],cmap="YlGnBu",
        square=True,
        linewidths=.5,vmax=400,
        xticklabels=3,#makes the labels every 3 nts
        cbar=True,cbar_ax=cbar_ax)
    ##plot the read length bar graph
    seaborn.barplot(x='index',y='total',
        data=df.reset_index(),ax=axs[2])
    ##now that we have plotted that data, we will compute the
    ##fraction in each frame.
    df[0]=df[0]/df['total']
    df[1]=df[1]/df['total']
    df[2]=df[2]/df['total']
    ##the next line performs several complicated things:
    ##df.reset_index() converts the index column to an proper column.
    ##the melt function reorganizes the data s.t. the frame columns
    ##are now rows. Doing these two things gets the data in a format
    ##for pandas to plot it.
    df2=pandas.melt(df.reset_index(),id_vars='index',
                    value_vars=[0,1,2],
                    var_name='Frame',
                    value_name='Frequency')
    ##now write the output
    seaborn.barplot(x='index',y='Frequency',hue='Frame',
                        data=df2,ax=axs[3])
    ##the next line moves the key outside of the figure
    axs[3].legend(bbox_to_anchor=(1.05, 1),title='Frame')
    ##add subplot axis title
    axs[0].set(xlabel='Position Relative Start Codon (nt)',
        ylabel='Read Length (nt)')
    axs[1].set(xlabel='Position Relative Stop Codon (nt)',
        ylabel='Read Length (nt)')
    axs[2].set(xlabel='Read Length',
        ylabel='Read Counts')
    axs[3].set(xlabel='Read Length')
    ##increase spacing between plots
    fig.tight_layout(pad=1.0)
    ##write the output--will make all three filetypes to allow the user
    ##to pick whichever they want.
    fig.savefig(f'{outPrefix}.png')
    fig.savefig(f'{outPrefix}.pdf')
    fig.savefig(f'{outPrefix}.svg')

def main(args):
    inFile,lowerBound,upperBound,upStart,downStart,upStop,downStop,\
        outPrefix=args[0:]
    ##first call the metaStartStop function
    dfStart,dfStop=metaStartStopHeatMap2.main([inFile,lowerBound,\
        upperBound,upStart,downStart,upStop,downStop,'blah'])
    ##now call the readLengthAndPhasingAnalysis2 function
    df=readLengthAndPhasingAnalysis2.main([inFile,lowerBound,upperBound,\
        'blah'])
    ##now plot everything
    mkPlot(dfStart,dfStop,df,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
