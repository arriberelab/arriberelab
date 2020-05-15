"""
Joshua Arribere, March 30, 2020

Script to make a heatmap of read length v number of reads.

Input: inFile.jam - .jam file
    readLengthLower, readLengthUpper - how short and long
        of read lengths you wish to examine
    startLeft, startRight, stopLeft, stopRight - how many
        nts to look up/downstream of start/stop. For ex:
        -10 50 -50 20 would be [-10,50], [-50,20]

Output: heatmap of reads per million

run as python metaStartStopHeatMap.py inFile.jam 15 18 
    -10 50 -50 20 outPrefix
"""
import sys, seaborn, matplotlib.pyplot
from logJosh import Tee
import pandas

def getPosition(txtList,index):
    """txtList is a list of txtInfo a la joshSAM. Will
    loop through all the txts. If the index position of
    the txts is all the same, then will return that
    position. Else will return 'na'. Will also require
    that all be Sense (not Antisense)."""
    transciptList = txtList.split('|')
    positions=[entry.split(':')[index] for entry in transciptList]
    positions=map(int,positions)
    positions=list(set(positions))
    if len(positions)==1:
        return positions[0]
    return 'na'

def parseJoshSAMToDataFrame(inFile,
        readLengthLower,readLengthUpper,
        leftBound,rightBound,startOrStop):
    """inFile is a .jam file. Will look between
    [leftBound,rightBound] relative to the startOrStop
    (1 if start, 2 if stop).
    Will initially add read counts for reads that
    are unambiguously assignable to a position on a txt.
    Then will normalize those to make rpm where m is million
    unambiguously assignable reads"""
    print('Restriction for uniquely mapping reads is on!')
    #for the next line you need to use data=0. Using
    #data=None initialized the df w/ NaN, which doesn't
    #work w/ subsequent incrementation.
    df=pandas.DataFrame(data=0,
        index=range(readLengthLower,readLengthUpper+1),
        columns=range(leftBound,rightBound+1))
    cntr=0
    with open(inFile,'r') as f:
        ##skip the first line
        next(f)
        ##now do analyses
        for line in f:
            if not line.startswith('@'):
                line=line.strip().split('\t')
                if len(line)>=10:
                    position=getPosition(line[9],startOrStop)
                    if position!='na' and line[8].endswith(':S'):
                        readLength=len(line[6])
                        if line[7]=='1:1' and line[3]=='-':
                            if readLength in range(readLengthLower,
                                                   readLengthUpper+1):
                                cntr+=1
                                if position in range(leftBound,
                                                     rightBound+1):
                                    df.loc[readLength,position]+=1
    #now normalize to make rpm
    norm=float(cntr)/1000000.
    for ii in range(readLengthLower,readLengthUpper+1):
        for jj in range(leftBound,rightBound+1):
            df.loc[ii,jj]/=norm
    #return the DataFrame
    return df

def mkHeatMaps(dfStart,dfStop,outPrefix):
    """Will plot the pair of dataframes dfStart and dfStop
    next to one another in outPrefix file"""
    #
    fig,axs=matplotlib.pyplot.subplots(nrows=2,
        ncols=1)
    ##subplot titles
    #axs[0].set_title('Position Relative Start Codon (nt)')
    #axs[1].set_title('Position Relative Stop Codon (nt)')
    #subplot of the start codon
    seaborn.heatmap(dfStart,ax=axs[0],cmap="YlGnBu",
        square=True,
        linewidths=.5,
        cbar_kws={"orientation": "horizontal",
                    "pad": 0.35,
                    "aspect": 40,
                    "label": 'RPM'})
    #subplot of the stop codon
    seaborn.heatmap(dfStop,ax=axs[1],cmap="YlGnBu",
        square=True,
        linewidths=.5,
        cbar_kws={"orientation": "horizontal",
                    "pad": 0.35,
                    "aspect": 40,
                    "label": 'RPM'})
    #
    #add subplot axis labels
    axs[0].set_xlabel('Position Relative Start Codon (nt)')
    axs[1].set_xlabel('Position Relative Stop Codon (nt)')
    axs[0].set_ylabel('Read Length (nt)')
    axs[1].set_ylabel('Read Length (nt)')
    #fig=axs.get_figure()
    #write output
    fig.savefig(f'{outPrefix}.png')
    fig.savefig(f'{outPrefix}.svg')

def main(args):
    #parse the in file
    inFile=args[0]
    #get the readLengths
    readLengthLower,readLengthUpper=map(int,args[1:3])
    #get the start/stop bounds
    startLeft,startRight,stopLeft,stopRight=map(
        int,args[3:7])
    if startLeft>0 or stopLeft>0:
        print("Head's up: if you want to look at positions \
            upstream of the start/stop codon, your \
            startLeft and stopLeft bounds need to be \
            negative. They currently are not.")
    #now that we have those inputs, it should
    #be possible to make the heatmap
    dfStart=parseJoshSAMToDataFrame(inFile,
        readLengthLower,readLengthUpper,
        startLeft,startRight,1)
    dfStop=parseJoshSAMToDataFrame(inFile,
        readLengthLower,readLengthUpper,
        stopLeft,stopRight,2)
    #and grab the output file name
    outPrefix=args[7]
    #make the heatmap
    #mkHeatMaps(dfStart,dfStop,outPrefix)
    ##return the df instead of plotting it
    return dfStart,dfStop

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
