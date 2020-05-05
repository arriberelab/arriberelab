"""
Joshua Arribere, April 6, 2020

Script to plot both read length and fraction in each frame for a library.

Input: file.jam - jam file
    lowerBound - smallest read size to consider
    upperBound - largest read size to consider

Output: two plots, one a histogram of readLength (x-axis) and read count
    (y-axis). Second plot will be readLength (x-axis) and fraction in
    frame (y-axis). Second plot will be a stacked bar graph.

run as python readLengthAndPhasingAnalysis.py inFile.jam X Y
    outPrefix
"""
import sys, common, pickle
from logJosh import Tee
import pandas
import seaborn,matplotlib

def getFrame(txtList,index):
    """txtList is a list of txtInfo a la jam. Will
    loop through all the txts. If the index position of
    the txts is all the same, then will return that
    position. Else will return 'na'. Will also require
    that all be Sense (not Antisense)."""
    SorAS=[entry.split(':')[3] for entry in txtList]
    SorAS=list(set(SorAS))
    if len(SorAS)==1 and SorAS[0]=='S':
        positions=[entry.split(':')[index] for entry in \
                txtList]
        positions=map(int,positions)
        positions=[entry%3 for entry in positions]
        positions=list(set(positions))
        if len(positions)==1:
            return positions[0]
    return 'na'

def getReadLengthAndFrameInfo(inFile,X,Y):
    """
    inFile is a jam file. Will consider reads in range [X,Y]
    inclusive. For each read in the jam file, will only focus on
    reads that are sense-stranded. Will build a DataFrame with indexes
    as readLengths and columns are 0, 1, or 2 frame. Will not tally
    reads--instead this information will be computed at the ends as a
    sum of the 0, 1, 2 frame read counts.
    """
    print('Restriction for uniquely mapping reads is on!')
    ##initialize the indexes
    #index=pandas.MultiIndex.from_tuples([0,1,2],
    #                                    names=['frame'])
    ##initialize the DataFrame
    df=pandas.DataFrame(data=0,
        index=range(X,Y+1),
        columns=range(3))
    ##loop through the file
    with open(inFile,'r') as f:
        ##ignore the first line
        next(f)
        ##now do something
        for line in f:
            if not line.startswith('@'):
                line=line.strip().split('\t')
                if line[7]=='1:1' and line[3]=='-':
                    frame=getFrame(line[9:],1)
                    if frame!='na':
                        readLength=len(line[6])
                        if readLength in range(X,Y+1):
                            df.loc[readLength,frame]+=1
    ##now compute the total
    df['total']=df[0]+df[1]+df[2]
    ##return
    return df

def mkPlots(df,outPrefix):
    """
    df is a pandas DataFrame where the first index column is an integer
    containing read length. The next three columns are titles 0, 1, 2 and
    are the frames of the reads. The last column is 'total' and is a sum
    of the 0, 1, 2 column read counts. Will plot the distribution of
    read lengths on the top, and the frequency of reads in each frame
    below.
    """
    ##first initialize the style
    seaborn.set_style("whitegrid", {'grid.linestyle': '--'})
    ##initialize matplotlib plot
    fig,axs=matplotlib.pyplot.subplots(nrows=2,ncols=1)
    ##plot the upper plot, of just read length v counts
    seaborn.barplot(x='index',y='total',
        data=df.reset_index(),ax=axs[0])
    axs[0].set(xlabel='Read Length',ylabel='Read Counts')
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
                        data=df2,ax=axs[1])
    axs[1].set(xlabel='Read Length')
    ##the next line moves the key outside of the figure
    axs[1].legend(bbox_to_anchor=(1.05, 1),title='Frame')
    ##increase spacing between plots
    fig.tight_layout(pad=1.0)
    ##save the figure
    fig.savefig(f'{outPrefix}.png')
    fig.savefig(f'{outPrefix}.svg')

def main(args):
    inFile,lowerBound,upperBound,outPrefix=args[0:]
    ##convert the bounds to integers
    lowerBound=int(lowerBound)
    upperBound=int(upperBound)
    
    ##extract the read length and fraction in frame info
    df=getReadLengthAndFrameInfo(inFile,lowerBound,upperBound)
    ##pickle the output to save it
    #with open(outPrefix+'.p','wb') as f:
    #    pickle.dump(df,f)
    
    ##unpickle the output to load it
    #with open(outPrefix+'.p','rb') as f:
    #    df=pickle.load(f)
    ##make the plots save
    #mkPlots(df,outPrefix)
    ##return the df instead of plotting it
    return df

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
