"""
Joshua Arribere, April 7 2020

Script to examine the distribution of Ribo-seq reads about
    positions in a group of transcripts.

Input: txtsAndPositions.csv - of the format:
    txtName,geneName,slipsite,ss_start\n
    where txtName and geneName are from ENSMEBL, slipsite is
    probably a 7mer, and ss_start is the start relative to the
    annotated start codon. There may be more columns than this, and
    a header will be required.
    genesTxtsCDSTxtLengths - of the format
        gene\ttxt\tCDSlength\ttxtLength
        Will be used to relate genes with transcripts and get
        CDS/txt lengths for normalization
    reads.joshSAM - Ribo-seq reads
    offsets.txt - of the format:
        readLength\toffset
        where a read of length readLength will have offset added to
        its 5'end position (relative to start codon)
    N - minimum number of reads/gene

Output: plot of read densities about the positions. May add other
    outputs as well (randomizations, phasing analysis, etc.).

run as python metaGeneByTxtPositions.py txtsAndPositions.txt
    genesTxtsCDSTxtLengths reads.joshSAM offsets.txt N outPrefix
"""
import sys, common, pickle
from logJosh import Tee
import pandas, numpy
import seaborn, matplotlib

def parsePositions(inFile):
    """
    Will parse a csv to a pandas DataFrame. Only requirement is there
    are columns named 'transcript' and 'gene', which will become the
    indexes.
    """
    with open(inFile,'r') as f:
        df=pandas.read_csv(f,
            dtype={'ss_start':"int64"},
            index_col=['transcript','gene'])
    #print(df),sys.exit()
    return df

def parseTxt(monoFile):
    """
    monoFile is of the format
    gene\ttxt\tCDSlength\ttxtLength
    Will parse to a DF
    """
    with open(monoFile,'r') as f:
        df=pandas.read_csv(f,
            sep='\t',
            index_col=['gene','txt'])
    return df

def parseOffsetDF(readOffsetsFile):
    """
    readOffsetsFile is of the format
    readLength\toffset
    Will parse to a DF
    """
    with open(readOffsetsFile,'r') as f:
        df=pandas.read_csv(f,
            sep='\t',
            index_col=['readLength'])
    return df

def passed(txtList,offset):
    txtNames=[entry.split(':')[0] for entry in txtList]
    pos1=[int(entry.split(':')[1])+offset for entry in txtList]
    pos2=[int(entry.split(':')[2])+offset for entry in txtList]
    strands=[entry.split(':')[3] for entry in txtList]
    if strands[0]=='S' and len(set(strands))==1:
        for entry1,entry2 in zip(pos1,pos2):
            if entry1<0 or entry2>-2:
                return 'na'
        return zip(txtNames,pos1,pos2)
    return 'na'

def parseReadFile(readsFile,offsetDF):
    """
    readsFile is a joshSAM file.
    offsetDF is a DF with readLength, offset
    This function will go through readsFile and create a dict of
    {txtID:{position:ct}} where ct has been normalized.
    """
    ##loop through reads
    aa={}
    with open(readsFile,'r') as f:
        for line in f:
            line=line.strip().split('\t')
            geneID=line[5]
            readLength=len(line[3])
            if readLength in offsetDF.index:
                ##
                offset=offsetDF.loc[readLength]['offset']
                ##
                numTxts=len(line[6:])
                txts=passed(line[6:],offset)
                ##
                if txts!='na':
                    if geneID not in aa:
                        aa[geneID]={}
                    for txtID,posRelStart,posRelStop in txts:
                        if txtID not in aa[geneID]:
                            aa[geneID][txtID]={}
                        if posRelStart not in aa[geneID][txtID]:
                            aa[geneID][txtID][posRelStart]=0
                        aa[geneID][txtID][posRelStart]+=1./numTxts
    return aa

def normalize(reads,txtDF,N):
    """
    reads is a dict of format:
    {gene:{txt:{position:ct}}}
    N is the read cutoff
    Will normalize the read counts for total reads per gene and
    CDS length.
    """
    print(f'Requiring at least {N} reads per gene')
    ##
    bb={}
    for gene in reads:
        for txt in reads[gene]:
            if (gene,txt) in txtDF.index:
                ##
                cdsLength=txtDF.loc[gene,txt]['CDSlength']
                ##
                totalReads=sum(reads[gene][txt].values())
                norm=totalReads/cdsLength
                if totalReads>=N:
                    ##
                    if gene not in bb:
                        bb[gene]={}
                    bb[gene][txt]=dict((k,v/norm) for k,v in \
                                reads[gene][txt].items())
    return bb

def computeMetaDFAboutPositions(reads,txtsAndPositions):
    """
    reads is a dict of format
    {gene:{txt:{position:ct}}}
    txtsAndPositions is a DataFrame
    Will return a DataFrame of position:ct
    where position is centered about the txtsAndPositions
    and ct is the average
    """
    ##this next line removes the periods from the txt and gene names
    txtsAndPositions.index=pandas.MultiIndex.\
        from_tuples([(x[0].split('.')[0],x[1].split('.')[0]) \
        for x in txtsAndPositions.index],
        names=['txt','gene'])
    ##sort indexes
    txtsAndPositions.sort_index(inplace=True)
    ##initialize the DataFrame
    ##the following number will be used to look this may nts up and
    ##downstream
    N=100
    df=pandas.DataFrame(data={'counts':0}, index=range(-N,N+1))
    ##
    cntr=0
    for index,row in txtsAndPositions.iterrows():
        txt,gene=index
        if gene in reads:
            if txt in reads[gene]:
                ss_start=row['ss_start']
                cntr+=1
                for ii in range(ss_start-N,ss_start+N+1):
                    if ii in reads[gene][txt]:
                        val=reads[gene][txt][ii]
                    else:
                        val=0
                    df.loc[ii-ss_start,'counts']+=val
    ##
    df=df.div(cntr)
    ##
    return df

def mkPlot(metaDF,outPrefix):
    """"""
    ##first convert the index to a column
    metaDF['index']=metaDF.index
    ##
    seaborn.set_style("whitegrid")
    ##
    ax=seaborn.lineplot(x='index',y='counts',data=metaDF)
    ##
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(3))
    ax.set_xlim(-30,30)
    ax.set(xlabel='Distance Relative SS Start',\
        ylabel='Normalized Ribo-seq Read Density')
    ##
    fig=ax.get_figure()
    fig.savefig(f'{outPrefix}.png')
    ##this next line is required to prevent that plot from appearing on
    ##top of another plot made later in this script
    matplotlib.pyplot.cla()

def reOrganize(metaDF):
    ##first slice things up
    zeroDF=metaDF.iloc[1::3,:]
    oneDF=metaDF.iloc[2::3,:]
    twoDF=metaDF.iloc[3::3,:]
    ##then reindex so they all have the same indexes now
    ##initially re-indexed as below, but then switched to make
    ##plotting easier later
    oneDF.index=range(-99,100,3)
    twoDF.index=range(-99,99,3)
    ##now make a DF with the total
    totalDF=pandas.DataFrame()
    totalDF['counts']=zeroDF['counts']+oneDF['counts']+twoDF['counts']
    totalDF['zero']=zeroDF['counts']/totalDF['counts']
    totalDF['one']=totalDF['zero']+oneDF['counts']/totalDF['counts']
    totalDF['two']=totalDF['one']+twoDF['counts']/totalDF['counts']
    ##subset the data to make plotting easier later b/c I cannot figure
    ##out how to reset xlim on the seaborn plot wihout it cutting the
    ##plot in two
    df2=totalDF.loc[-30:30]
    ##
    return df2

def mkBarChartStacked(metaDF,outPrefix):
    """
    """
    ##first reorgnize the metagene to display frame
    metaDF=reOrganize(metaDF)
    metaDF['index']=metaDF.index
    ##
    ax=seaborn.barplot(x='index',y='two',data=metaDF,
        color=(0,114/255,178/255))
    ax.set_xlim(-10,10)
    seaborn.barplot(x='index',y='one',data=metaDF,
        color=(0,158/255,115/255))
    seaborn.barplot(x='index',y='zero',data=metaDF,
        color=(240/255,228/255,66/255))
    ##
    ax.set_ylim(0,1)
    ax.set(xlabel='Distance Relative SS Start',\
        ylabel='Fraction of Normalized Ribo-seq Read Density in Frame')
    ##
    fig=ax.get_figure()
    fig.savefig(f'{outPrefix}.bar.png')

def main(args):
    ##first parse the input
    txtsAndPositionsFile,txtGeneFile,\
        readsFile,readOffsetsFile,N,outPrefix=\
        args[0:]
    ##
    outPrefix+='_'+N
    N=int(N)
    ##parse the monoTxt file
    txtDF=parseTxt(txtGeneFile)
    ##parse the offset file
    offsetDF=parseOffsetDF(readOffsetsFile)
    ##parse the reads
    reads=parseReadFile(readsFile,offsetDF)
    ##pickle the reads
    with open(f'{outPrefix}.p','wb') as f:
        pickle.dump(reads,f)
    ##unpickle the reads
    with open(f'{outPrefix}.p','rb') as f:
        reads=pickle.load(f)
    ##normalize the read counts
    reads=normalize(reads,txtDF,N)
    ##parse the file of txts and positions
    txtsAndPositions=parsePositions(txtsAndPositionsFile)
    ##now get the metaGene about those positions
    metaDF=computeMetaDFAboutPositions(reads,txtsAndPositions)
    ##plot the output
    mkPlot(metaDF,outPrefix)
    ##pickle the reads
    with open(f'{outPrefix}.2.p','wb') as f:
        pickle.dump(metaDF,f)
    ##unpickle the reads
    with open(f'{outPrefix}.2.p','rb') as f:
        metaDF=pickle.load(f)
    ##make a stacked bar chart of the output
    mkBarChartStacked(metaDF,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
