#!/usr/bin/env python3
# Name: Nitin Vidyasagar (nvidyasa)


"""
200704, Nitin Vidyasagar

Represents fractional read abundances for genes for Start and Stop Codon

Input: jamFilesList.txt, favGenes.csv
Output: 

jamFilesList.txt is a tab delineated text file
inFile1\tdirectory

favGenes.csv is formatted as follows:

Run as python3 fractionalReadAbundance.py files.text favGenes.csv outPrefix




"""
import os
import sys, csv, time
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
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
        self.parser = argparse.ArgumentParser(description="*.jam files from assignReadsToGenes5.py")
        self.parser.add_argument("inFile1", action="store", nargs="+", help=" tab delineated text file for jam files")
        self.parser.add_argument("inFile2", action="store", nargs="+", help="csvFile for genes")
        self.parser.add_argument("outPrefix", action="store", nargs="+",
        help="Outprefix for pdf for Graphs")
   

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class HistoneContainers():
#     """ Creates an object that has lists for the histone genes """
    def __init__(self):
#         """Initializes lists for the histone genes and their dictionaries for the counts"""
        self.H2AList=[]
        self.H2BList=[]
        self.H3List=[]
        self.roughH4List=[]
        self.finalH4List=[]
        
        #Initialize the dictionaries as attributes for this class
       
        

    def csvParser(self, csvFile):  
#         """ appends columns from csv to each of the histone gene lists"""
        with open(csvFile,'r') as csv_file: #https://stackoverflow.com/questions/48692264/how-to-import-csv-column-as-list-in-python
            lines = csv_file.readlines()
            for line in lines:
                data = line.split(',')
                self.H2AList.append(data[0])
                self.H2BList.append(data[1])
                self.H3List.append(data[2])
                self.roughH4List.append(data[3])
                
            # Removes the first item of the list which is just H2A, H2B, H3. Remove last
            # item too of H4 respectively
            self.H2AList.pop(0)
            self.H2BList.pop(0)
            self.H3List.pop(0)
            self.roughH4List.pop(0) #strangely, got new line \n characters. Need to strip them out
            self.roughH4List.pop()
            
            for geneID in self.roughH4List:   # Strips newLine characters from this list and appends to the final H4 list
                self.finalH4List.append(geneID.strip()) 

        return self.H2AList, self.H2BList, self.H3List, self.finalH4List


    def histoneDictionaryCreator(self):
        """Method appends WB Gene name from the histone genes as keys to the dictionary. Sets initial values of the keys as 0"""

        
        self.H2ACounts= {gene: 0 for gene in self.H2AList}
        self.H2BCounts= {gene: 0 for gene in self.H2BList}
        self.H3Counts= {gene: 0 for gene in self.H3List}
        self.H4Counts= {gene: 0 for gene in self.finalH4List}

        return self.H2ACounts, self.H2BCounts, self.H3Counts, self.H4Counts
            
            
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

 
        




    def countsOfGenes(self,fileName1,H2AList,H2BList,H3List,H4List):
        """
        Collects counts of sense reads of genes in the inFile.jam file and also the metagene counts/fraction of reads over start codon/stop codon

        Input:  None
        Output: updates dictionary of the genes with the counts and the total geneCounts

        This function will weight genes that map to multiple areas:
        eg. Count of a Gene that mapped to 5 areas=  1/5 for each time the read was mapped
        Count of gene that was uniquley mapped is = 1/1 

        To do this, we will first for each line identify what the gene the read mapped to was. We create keys in the dictionary
        that are the WBgenes. 

        We will only be counting sense reads in our gene counts for each gene. However, we will also count Anti-sense reads when we consider
        the total reads of the inFile.jam

        To deal with multiply mapping genes for the metagene counts, I decided to just add them into the atgReads and collapse them at this
        step rather than in the downstream steps. By collapse, I mean combine the counts into one gene.

        For ATG and STOP reads, we want to make sure that all the transcript info have the same range of position from the ATG and STOP codon

        """
        atgReads={'H2A':0, 'H2B':0, 'H3':0,'H4':0,'EEF-2':0,'UNC-54(NSD)':0}
        stopReads={'H2A':0, 'H2B':0, 'H3':0,'H4':0,'EEF-2':0,'UNC-54(NSD)':0}
        geneCounts={'H2A':0, 'H2B':0, 'H3':0,'H4':0,'EEF-2':0,'UNC-54(NSD)':0}
        totalReadCounts=0
   




        with open (fileName1, 'r') as f:
            for line in f:
                linesSplit=line.split() #creates a list of all elements in the file
                multipleMapData= linesSplit[5]  #access number of times read was mapped in genome. Format is x:y where y represents the number of different genes  read was mapped to
                wbGeneID= linesSplit[6] #access the wbGeneID
                lengthRead=int(linesSplit[4]) #length of read
                listForMappedTimes=multipleMapData.split(":") # split string to acces 'x' and 'y' using indices
                numberOfDifferentGenesMapped= int(listForMappedTimes[1]) # the we have accessed and stored y which is the genes the single read mapped to.
                
                #we now want to check and ensure that we are accessing sense reads only: ie. the transcripts the read mapped to are all S
                transcriptList= linesSplit[7:]  
                senseAntisenseList=[]


                for transcript in transcriptList:
                    elementsOfTransplit=transcript.split(":")
                    senseAntisense=elementsOfTransplit[3]
                    senseAntisenseList.append(senseAntisense)


                if "AS" not in str(senseAntisenseList):

                    if lengthRead in range(15,19):
                    	#We are calculating the total number of reads for each of the histone metagenes now. We will count only 15-18 nt reads for the metagene
                        # if wbGeneID == "WBGene00001167": #This is EEF-2 a gene not affected by pelo and ski
                        if any("F25H5.4" in transcript for transcript in transcriptList):#I want to know if the gene can be multiply mapping
                            geneCounts['EEF-2']+=float(1/numberOfDifferentGenesMapped)
                        if any("unc-54" in transcript for transcript in transcriptList):
                            geneCounts['UNC-54(NSD)']+=float(1/numberOfDifferentGenesMapped)
                        elif wbGeneID in H2AList:
                            geneCounts['H2A']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes
                        elif wbGeneID in H2BList:
                            geneCounts['H2B']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes
                        elif wbGeneID in H3List:
                            geneCounts['H3']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different gene
                        elif wbGeneID in H4List:
                            geneCounts['H4']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes

                        ATGList=[] #We will append all the start ints to this list
                        STOPList=[]#We will append all the stop ints to this list

                        for transcript in transcriptList:
                            elementsOfTransplit=transcript.split(":")
                            ATG=int(elementsOfTransplit[1]) #int value of the distance from start codon
                            ATGList.append(ATG)

                            STOP=int(elementsOfTransplit[2])
                            STOPList.append(STOP)

                        if all(-12<x <40 for x in ATGList):# if all the transcripts have reads less than 30
                            # if wbGeneID == "WBGene00001167":
                            if any("F25H5.4" in transcript for transcript in transcriptList):
                                atgReads['EEF-2']+=float(1/numberOfDifferentGenesMapped)
                            if any("unc-54" in transcript for transcript in transcriptList):
                                atgReads['UNC-54(NSD)']+=float(1/numberOfDifferentGenesMapped)
                            elif wbGeneID in H2AList:
                                atgReads['H2A']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes
                            elif wbGeneID in H2BList:
                                atgReads['H2B']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different gen
                            elif wbGeneID in H3List:
                                atgReads['H3']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different gen
                            elif wbGeneID in H4List:
                                atgReads['H4']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes

                        if all(-12<x < 40 for x in STOPList):# if all the transcripts have reads less than 30
                            # if wbGeneID == "WBGene00001167":
                            if any("F25H5.4" in transcript for transcript in transcriptList):
                                stopReads['EEF-2']+=float(1/numberOfDifferentGenesMapped)
                            if any('unc-54'in transcript for transcript in transcriptList):
                                stopReads['UNC-54(NSD)']+=float(1/numberOfDifferentGenesMapped)

                            elif wbGeneID in H2AList:
                                stopReads['H2A']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes
                            elif wbGeneID in H2BList:
                                stopReads['H2B']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes
                            elif wbGeneID in H3List:
                                stopReads['H3']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes
                            elif wbGeneID in H4List:
                                stopReads['H4']+=float(1/numberOfDifferentGenesMapped)#we add a fractional read count if the read maps to multiple different genes








                #For funsies, we calculate all the reads present in the file just incase anyone decides to do RPM normalization

                if numberOfDifferentGenesMapped==1:
                    totalReadCounts+=1
                else:
                    totalReadCounts+= float(1/numberOfDifferentGenesMapped)

                
               
            	
            return totalReadCounts, geneCounts, atgReads,stopReads

class Seaborn:

	def barGraph(self,dataframe1,dataframe2,outPrefix):
		"""Creates a bargraph using seaborn"""
		# Set up the matplotlib figure
		
		sns.set(style="white")
		plt.ylim(0, 1)

		ax = sns.barplot(x="Library", y="FractionReads", hue="Protein", data=dataframe1)
		ax.set_xlabel('Genotype', fontsize=16)
		ax.set_ylabel('Fraction Counts', fontsize=16)
		ax.set_title("Fraction Over START", fontsize=16)
		plt.savefig(outPrefix+"FractionStart.pdf", dpi=300)
		plt.close()


		sns.set(style="white")
		plt.ylim(0, 1)
		ax2 = sns.barplot(x="Library", y="FractionReads", hue="Protein", data=dataframe2)
		ax2.set_title("Fraction Over STOP",fontsize=16)
		ax2.set_xlabel('Genotype', fontsize=16)
		ax2.set_ylabel('Fraction Counts', fontsize=16)
		plt.savefig(outPrefix+"FractionSTOP.pdf", dpi=300)
		plt.close()

		













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

    jamParse= JoshSamParser()#create object from JoshSamParser class from the command line
    histones= HistoneContainers()
    readFiles=myCmdLine.args.inFile1[0]
    csvFile=myCmdLine.args.inFile2[0]
    outPrefix=myCmdLine.args.outPrefix[0]

    sb= Seaborn()

    #generate the histone lists
    H2AList,H2BList,H3List,H4List= histones.csvParser(csvFile)

    

    pd.set_option("display.max_rows", None, "display.max_columns", None) #allows me to remove any limits on number of rows and columns to view
    
    
   
    

 



    ATGList=[]
    STOPList=[]

    if os.path.isfile(readFiles):
        with open(readFiles,'r') as f:
            freader=csv.reader(f,delimiter='\t')
            for row in freader:
                name,fileName=row[0],row[1]
       
                tmpFileName=fileName
                print ('Working on ' + name)
                totalReadCounts, geneCounts,atgReads, stopReads= jamParse.countsOfGenes(tmpFileName,H2AList,H2BList,H3List,H4List)

                totalReadCountsinMillions= totalReadCounts* 10**(-6)

                #Now we will calculate the fraction of reads at ATG and TAA by dividing by geneCounts for each wbGeneID
                print ('Working on fraction ATG reads')
                for metagene, value in atgReads.items():
                	atgReads[metagene]=atgReads[metagene]/geneCounts[metagene]

                print ('Working on fraction Stop reads')
                for metagene, value in stopReads.items():
                	stopReads[metagene]=stopReads[metagene]/geneCounts[metagene]


                print ('Converting dictionary to dataframe and storing in ATG list')
                atgtempdf=pd.DataFrame.from_dict(atgReads, orient='index')
                atgtempdf=atgtempdf.reset_index()
                atgtempdf=atgtempdf.rename(columns={0:"FractionReads","index":"Protein"})
                atgtempdf['Library']=name
                ATGList.append(atgtempdf)
              
                print ('Converting dictionary to dataframe and storing in STOP list')
                tempdf=pd.DataFrame.from_dict(stopReads, orient='index')
                tempdf=tempdf.reset_index()
                tempdf=tempdf.rename(columns={0:"FractionReads","index":"Protein"})
                tempdf['Library']=name
                STOPList.append(tempdf)

    print ('Concatenating all dataframes together for STOP')
    STOPdf=pd.concat(STOPList)
    print(STOPdf)
    print ('Concatenating all dataframes together for START')
    ATGdf=pd.concat(ATGList)
    print(ATGdf)

    print ('Dataframes are ready! Let us Graph')
    sb.barGraph(ATGdf,STOPdf,outPrefix)
    print ('The graphs are ready! Have a wonderful day, ^_^')


if __name__ == '__main__':
    Tee()
    main()