#!/usr/bin/env python3
# Name: Nitin Vidyasagar (nvidyasa)
# Group Members: Hannah Knotter

"""
200409, Nitin Vidyasagar

Program purpose: To loop through a Joshjam file and return a pkl object for the dataframe of all genes and their RPM normalized counts. Only takes count of sense reads

Input: inFile.jam  (format from wrapper9)
Output: pickled dataframe for gene counts 

Run as python3 multipleMappedGeneCounts.py inFile.jam outFile.pkl

Required modules:
  
sys, csv, time
pandas as pd
pickle
logJosh import Tee

"""
import numpy as np
import sys, csv, time
import pandas as pd
import pickle
from logJosh import Tee






class CommandLine():
    """Handle the command line to determine inFiles and outFile.
    Only arguments that are specified in the command line are the file names."""
    def __init__(self, inOpts=None):
        """Implement a parser to interpret command line argv string using argparse."""
        import argparse
        self.parser = argparse.ArgumentParser(description = 'parses *.jam files from assignReadsToGenes5.py. Since we use RPM normalization, you should use unfiltered reads ', 
                                             
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = 'python3 %(prog)s <input >output'
                                             )
        self.parser = argparse.ArgumentParser(description="*.redudantAndUnique.jam files from wrapper9")
        self.parser.add_argument("inFile1", action="store", nargs="+", help="*.jam generated from assignReadsToGenes5.py")
        self.parser.add_argument("outFile1", action="store", nargs="+",
        help="output file is pickled data frame for genes counts: outFile1.pkl")
   

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

            
            
class JoshSamParser:
    """
    Creates an object to parse joshSamFile and contains methods to collect counts of genes.

    initialized attributes:
    fileName1: this is the jam file. The reads are mapped to genes. If read maps to multiple parts of genome, it is given 
    a count of 1/n where n is the number of different genes the read mapped to
    totalReadCounts  
    
    methods:
    countsOfGenes: This program caculates the counts of genes in joshJam file



    """

    def __init__ (self,fileName1):
        self.fileName1= fileName1
        self.totalReadCounts=0
        self.geneCountsDict={}



    def countsOfGenes(self):
        """
        Collects of sense reads of genes in the inFile.jam file

        Input:  None
        Output: updates dictionary of the genes with the counts and the total geneCounts

        This function will weight genes that map to multiple areas:
        eg. Count of a Gene that mapped to 5 areas=  1/5 for each time the read was mapped
        Count of gene that was uniquley mapped is = 1/1 

        To do this, we will first for each line identify what the gene the read mapped to was. We create keys in the dictionary
        that are the WBgenes. 

        We will only be counting sense reads in our gene counts for each gene. However, we will also count Anti-sense reads when we consider
        the total reads of the inFile.jam

        """


        with open (self.fileName1, 'r') as f:
            next(f) #skip first line as this is header line from wrapper9 
            for line in f:
                linesSplit=line.split() #creates a list of all elements in the file
                multipleMapData= linesSplit[7]  #access number of times read was mapped in genome. Format is x:y where y represents the number of different genes  read was mapped to
                geneInfo= linesSplit[8].split(":") #access the wbGeneID
                geneID= geneInfo[0]
                senseAntisense= geneInfo[1]
                lengthRead=len(linesSplit[6]) #length of read
                listForMappedTimes=multipleMapData.split(":") # split string to acces 'x' and 'y' using indices
                numberOfDifferentGenesMapped= int(listForMappedTimes[1]) # the we have accessed and stored y which is the genes the single read mapped to.

                #we now want to check and ensure that we are accessing sense reads only: ie. the transcripts the read mapped to are all S


                if senseAntisense=="S":

                    if lengthRead in range(15,19):
                    	if geneID not in self.geneCountsDict.keys():
                    		self.geneCountsDict.update({geneID:float(1/numberOfDifferentGenesMapped)}) #we add a fractional read count if the read maps to multiple different genes
                    	elif geneID in self.geneCountsDict.keys():
                    		self.geneCountsDict[geneID]+=float(1/numberOfDifferentGenesMapped)


                if numberOfDifferentGenesMapped==1:
                    self.totalReadCounts+=1
                else:
                    self.totalReadCounts+= float(1/numberOfDifferentGenesMapped)

                # my logic is to append a wbGENE as a key if not in the dictionary and set value as the 1/n where n is the number of different mapped regions
               
            	
            return self.geneCountsDict, self.totalReadCounts
               
       

def main(inCL=None):
    """
    defines the inputs to the program's command line. Creates dataframe of all histone genes and their counts
    Input: arguments from command line
    Output: printed dataframes and pickled dataframes
    """

    if inCL is None:
        myCmdLine = CommandLine()
    else:
        myCmdLine = CommandLine(inCL)


    jamFile= JoshSamParser(myCmdLine.args.inFile1[0]) #create object from JoshSamParser class from the command line
    
    
    multipleMappedGeneCounts, totalReadCounts= jamFile.countsOfGenes()

    #convert total read counts to Millions

    totalReadCountsinMillions= totalReadCounts* 10**(-6)
    
    #I will be using these functions to create a dataframe from dictionary. The keys (WBgene or the meta-gene of H2A,H2B,H3 and H4) will be index but the values (counts) will be columns but I am putting library name
    geneCountDataFrame= pd.DataFrame.from_dict(multipleMappedGeneCounts, orient='index').rename(columns={0:myCmdLine.args.inFile1[0]})

    #This step calculates the RPM normalized values for a library
    geneCountDataFrame[myCmdLine.args.inFile1[0]]= geneCountDataFrame[myCmdLine.args.inFile1[0]].div(totalReadCountsinMillions)

   

    pd.set_option("display.max_rows", None, "display.max_columns", None) #allows me to remove any limits on number of rows and columns to view

    


    print(geneCountDataFrame)


    return geneCountDataFrame.to_pickle(myCmdLine.args.outFile1[0])
    
if __name__ == '__main__':
    Tee()
    main()