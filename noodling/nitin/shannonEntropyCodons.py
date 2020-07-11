"""
Nitin Vidyasagar,June 30th, 2020

Script analyzes Coding Sequences aligned with "-" character fasta format to determine the variability of codons at that position

Input: fasta file of aligned CDS sequences (I would use https://mafft.cbrc.jp/alignment/software/macportable.html)
Output: jpg of graph and *.txt of highly conserved codons
"""

import sequenceAnalysis as seqan 
import csv
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import seaborn as sns


class CommandLine():
    """Handle the command line to determine inFiles and outFile.
    Only arguments that are specified in the command line are the file names."""
    def __init__(self, inOpts=None):
        """Implement a parser to interpret command line argv string using argparse."""
        import argparse
        self.parser = argparse.ArgumentParser(description = 'aligned fasta files to plot shannon entropy line graph'
                                             
                                             ,add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = 'python3 %(prog)s <input >output'
                                             )
        self.parser = argparse.ArgumentParser(description="*.redudantAndUnique.jam files from wrapper9")
        self.parser.add_argument("inFile1", action="store", nargs="+", help="*.aligned.fa from mafft")
        self.parser.add_argument("outPrefix", action="store", nargs="+",
        help="saved image for graph")
   

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)



def main(inCL=None):
    ''' 
    '''
    if inCL is None:
        myCmdLine = CommandLine()
    else:
        myCmdLine = CommandLine(inCL)

    fasta=seqan.FastAreader(myCmdLine.args.inFile1[0])# create fasta object
    nuc= seqan.NucParams() #create unique object
    outPrefix= myCmdLine.args.outPrefix[0]

    conservedText=outPrefix + str('.MinAndMaxConservationCodons.txt')
   
    #basicallly, if codon position 
    
    for head,seq in fasta.readFasta(): #iterate through fasta file 
        
        codonPos=nuc.codonTypesList(seq) #return the codon number and aa dictionary for sequence. Assume sequences are same length
    
    #We now have a dictionary where we have codon positon as key and all the possible values that codon one could occupy: eg {1: [CODON1,CODON1,CODON1]}_
    
    probabilityDict=nuc.probabilityDictCreator()
    

    # Probability is defined as follows: Lets say codon 1 across different sequences is ATG  for 4 sequences and TTG for 3 sequences
    # probabilities set at codon one is Number of ATG/7


    #We will now calculate the Shannon Entropy
    shannonEntropyDict=nuc.shannonEntropyCalculator(probabilityDict)


    print(shannonEntropyDict)
    minShannonEntropy=min(shannonEntropyDict.values())
    maxShannonEntropy=max(shannonEntropyDict.values())

    

    
    

    

    with open(conservedText,'w') as f:
        
        f.writelines('Here are the most conserved codons positions and their codons + probabilities\n')
        

        for key, value in shannonEntropyDict.items():
            if shannonEntropyDict[key]==minShannonEntropy:
                f.writelines('Codon '+str(key)+' had these codons ' + str(probabilityDict[key])+"\n")

        f.writelines('\n')

        f.writelines('Here are the least conserved codons and their positions + probabilities\n')

        for key, value in shannonEntropyDict.items():
            if shannonEntropyDict[key]== maxShannonEntropy:
                f.writelines('Codon '+str(key)+' had these codons ' + str(probabilityDict[key])+"\n")



                
            


    #create dataframe
    df=pd.DataFrame.from_dict(shannonEntropyDict,orient='index')
    df.reset_index(inplace=True)
    df=df.rename(columns = {"index":"CodonPosition",0:'ShannonEntropy'})
    print(df)


    # #Now we plot the data
    sns.set()
    plt.figure(figsize=(5,5))#we set the size of the figure
    ax = sns.lineplot(x="CodonPosition", y="ShannonEntropy", data=df)
    plt.savefig(outPrefix + ".JPEG")
    








    

   















main()

	





