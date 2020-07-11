"""
#!/usr/bin/env python3
# Name: Nitin Vidyasagar



The script can analyze a DNA/RNA sequence and Protein sequence from a fasta input
and compute statistics on the sequences. This program also contains classes to 
calculate details on Open reading frames

Classes include 


1) OrfFinder
This class contains methods to generate a reverse complement and find information
(frame, positionStart, positionStop and Length)
of open reading frames from a sequence file of dna sequence

2)NucParams 
This class contains methods to calculate
the nucleotide counts, codon composition, parse strings to remove 
unwanted charactes

3) ProteinParam
This class contains methods to calculate
the pI, molecular weight, composition and 
molar extinction, mass extinction of an amino acid sequence 
from an input string of protein sequence


4)FastAReader
This class contains methods to parse a fasta File and retrieve
header and sequence informations. 



Modules required:
sys



"""





import sys
from collections import defaultdict
import math

class OrfFinder:
	"""
	This class contains methods to generate a reverse comelement of DNA sequence and 
	collect information on ORFs such as the reading frame, posiitonofStart, PositionofStop and length
    
    Input: startCodonList from commandline, stopcodonlist from commandline and longestGene condition
    Output: Object that can compute ORF details and generate reverse complement


    class attributes:
    self.geneList=[] #This list will house the master information about gene candidates
	self.startCodon=set(startCodon)  #converts input list of starts to consider as a set
	self.stopCodons=set(stopCodon) #converts input list of stops to consider as a set
	self.longestGene=longestGene # This is a parameter that is default of false unless specified which is true. 
    

    To instantiate: myOrf= OrfFinder(startCodon,Stopcodons, longestGene)

    initialized: list for ORF gene candidates, set of start codons, set of stop codons and boolean of longestgene (true or false depending on
    commandline argyment)

    methods: geneFinder(dnaSequence, reverseComplement)

    Required module: sys



	"""

	def __init__(self,startCodon,stopCodon,longestGene):
		""" 
		Initializes lists for orfcandidates, set of start codons and stop codons to consider and the boolean for longestGeme and 
		a dictionary for generating a reverse complement. 
		
		input: start codon List, Stop codon list and longest gene, all from command line arguments
		output: Initialized lists for orfcandidates, set of start codons and stop codons to consider and the boolean for longestGeme

		"""
		self.geneList=[] #This list will house the master information about gene candidates (list of lists)
		self.startCodon=set(startCodon) #create a set of codons to consider as starts
		self.stopCodons=set(stopCodon) #create a set of codons to consider as starts
		self.longestGene=longestGene #this evaluates to true if specified on the command line (from findORFs)
		self.complementDict= {'A':'t', 'G':'c','U':'a','T':'a','C':'g'} # I will use this dictionary to generate a reverse complement


		



	def geneFinder(self,dnaSequence, reverseComplement):
		"""
		Scans a string of DNA sequence and its reverse complement and appends open reading frame information (length, reading frame and locations of start and stop in top strand coordinates)
		as a list in the self.geneList attribute.

		Input: dnasequence and reverseComplement (both strings from fasta files and we generated reverse complement using method of OrfFinder Class
		Output: geneList updated with lists containing information of ORFS


		"""

		if dnaSequence:	 #If we consider the top strand, we need to apply a different set of calculations on it than the reverse complement
			
			for frame in range(0,3): #for loop for every frame
				ts=[0] #this is a running list of starts for each frame. By setting ts=[0] we assume there is a start upstream of given sequence, hence helpful for dangling stop condition
				
				
				for p in range(frame,len(dnaSequence),3):
					
					codon= dnaSequence[p:p+3] #we created codons (groups of 3 letters) from the reading frame as starting position
			
					
					if codon in self.startCodon: #we not populate our list of start codons
						if p==frame:   #clear the list if there is a start codon in first position of each frame
							ts=[]
						ts.append(p) 


					elif codon in self.stopCodons:
						#This will generate gene information for all types of starts
						#clear list whenever we find a stop codon 
						for start in ts:
							positionStop=p+3
							positionStart= start + 1
							length= positionStop - positionStart +1
							gene= [(frame+1), positionStart, positionStop, length]  #the format for recording a potential gene of top strand is [positionStart, positionStop, lengthGene, readingFrame]
							self.geneList.append(gene)
							if self.longestGene==True:
								break #to stop the foreloop from going past first start position in ts because the first element is the longest gene
						ts=[] #clear list of starts each time we encounter a stop so we do not have open reading frames with stop codons in them
						

						

			# USE THIS FOR NO START OR STOP IN FRAME and For Dangling starts. We enter this part of the code once we hit the last stop and it was tke care of.
				for start in ts:
					positionStart= start+1
					positionStop= len(dnaSequence) # we assume end of dna sequence is the stop position
					length= positionStop-positionStart+1
					gene= [(frame+1), positionStart, positionStop, length]
					self.geneList.append(gene)
					if self.longestGene == True:
						break
				ts=[]

			

		if reverseComplement:
			
			for frame in range(0,3): #for loop for every frame
				ts=[0]
				for p in range(frame,len(reverseComplement),3):
					codon= reverseComplement[p:p+3]
					if codon in self.startCodon:
						if p==frame:   #clear the list if there is a start codon in first position of each frame
							ts.pop(0)
						ts.append(p)

					elif codon in self.stopCodons:
						for start in ts:
							positionStop=(len(reverseComplement)-(start)) #the "stop" position is the end of the gene on the top strand coordinates. The corresponds to the location of start codon on top strand
							positionStart= (len(reverseComplement)-(p+2)) #this is the position of the start codon 
							length= positionStop - positionStart +1
							gene= [-(frame+1), positionStart, positionStop, length]  #the format for recording a potential gene of top strand is [positionStart, positionStop, lengthGene, readingFrame]
							self.geneList.append(gene)
							if self.longestGene is True: #since the longest gene for a stop codon is the first start in ts list, we break the algorithm after first iteration
								break
						ts=[] #clear the start lists so we don't have open reading frames with stops in them.
				

				# USE THIS FOR NO START OR STOP IN FRAME and For Dangling starts. We enter this part of the code once we hit the last stop and it was tke care of.
				for start in ts:
					
					positionStop= (len(reverseComplement)-start)  #the location of the end of the ORF (hence, position stop) is the location of first letter of start codon on top strand
					positionStart= 1  #for dangling starts, this is always going to be the start position on the top strand
					length= positionStop-positionStart+1
					gene= [-(frame+1), positionStart, positionStop, length]
					self.geneList.append(gene)
					if self.longestGene is True:#since the longest gene for a stop codon is the first start in ts list, we break the algorithm after first iteration
						break
				ts=[]
		return self.geneList
						

	def reverseComplementCreator(self,dnaSequence):
		"""
		Creates the reverse complement by converting nucleotides of an input string of DNA with values from self.complementDict. Upper case and flip string to 
		create the reverse complement 
		
		Input: DNA Sequence from fasta file
		Output: reverecomplement of input dna string

		"""
		
		
		for orginalNuc in dnaSequence: #convert every nucleotide to complement
		    if orginalNuc in self.complementDict.keys(): #check if nucleotide is key in the complementDict
		        dnaSequence=dnaSequence.replace(orginalNuc,  self.complementDict[orginalNuc])  #if it is, replace with lower case complement
		        
		  
		reverse=dnaSequence[::-1] #flip the string (reverse)
		self.reverseComplement=reverse.upper() #convert to upper case
		return self.reverseComplement

class NucParams:
    """
    This class contains methods to compute 
    the nucleotide counts, codon composition, addSequence
    
    Input: String of amino acid sequence
    Output: Object that can compute pI, molecular weight, composition and 
    molar extinction, mass extinction of an amino acid sequence
    
    class attributes:
    self.rnaCodonTable
    self.dnaCodonTable
    self.nucComp
    self.aaComp
    self.nucCountValue

    To instantiate: myNuc= NucParams(sequence)

    initialized: dictionaries for aaComposition, nucleotide compotion,
            codon composition, and a nucCountValue integer

    methods:  addSequence(inSeq), aaComposition(cleanedSeq), nucComposition(cleanedSeq),
    codonComposition(cleanedSeq), codonListGenerator(cleanedSeq), nucCount(cleaned_sequence)

    Required module: sys


    
    """
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    
       
    

    def __init__ (self, inString=''):
        """
        Creating dictionaries for methods to append and utilize

        input: string of DNA Sequence/RNA sequence from fasta File
        output: initalized dictionaries for codon composition,
        nucleotide composition, amino acid composition and nucleotide
        count value.

        """
        
        self.nucComp={'A':0, 'T':0 , 'G':0, 'C':0, 'U':0, 'N':0} #sets values of valid nucleotides to 0
        self.codonComp= {key:0 for key,value in self.rnaCodonTable.items()} #sets values of valid codons to 0
        self.aaComp = {
        'A': 0,  'G': 0,  'M': 0, 'S': 0, 'C': 0,
        'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
        'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
        'W': 0,  'F': 0, 'L': 0, 'R': 0, 'Y': 0,
        '-':0
        }
        self.nucCountValue=0 
        self.codonListofLists=[]
        self.counter=0 #this will allow for a special set of operations when I run codon
        self.probabilityDict={}
        

        
    def addSequence (self,inSeq):
        """
        Parses a string sequence of DNA/RNA input to return uppercase letters and no whitespace. Also eliminates 
        unwanted characters from the sequence string. Returns a string containing
        only valid letters {A, T, G, U, N} in upper case.

        input: sequence string from parsed fasta File
        output: cleaned up, uppercase string of nucleic acid sequence
        """
        upperSequence =inSeq.upper() #convert string to upper case

        upperSequence.replace(" ","") #get rid of white space


        emptyPlaceHolder='' #where we will add the edited sequence to

        #Checking if letters are valid nucleotides and adding only valid nucleotides to the emptyPlaceHolder string
        for index in range (0, len(upperSequence)): #for loop to access all indices of string character
            if upperSequence[index] in self.nucComp.keys(): #if character is present as a valid key in self.nucComp Dict
                emptyPlaceHolder+= str(upperSequence[index]) #add the character to the emptyPlaceHolder String

            else:
                emptyPlaceHolder+= "" 

        self.cleanedSequence= emptyPlaceHolder #finalized string for us to use
        return self.cleanedSequence
        
       
            
            
    def aaComposition(self):
        """ 
        Returns aaComp Dictionary with counts of each amino acid.
        Uses output from codonListGenerator
        to get the list of codons. Then, the aaComp dictionary is returned by checking
        if codons in sequence are keys in RNA Codon Table. If codon is prsent as a key (ie. valid codon)
        update the value of the amino acid  key in aaComp coded by that codon

        Input: String of cleaned DNA Sequence
        Output: Updated aaComp dictionary with counts

        """
        
        codonList= self.codonListGenerator() #generate list of codons from sequence
        for codon in codonList: #for loop to check is codon is present in the rnaCodonTable keys
            if codon in self.rnaCodonTable.keys():
                aminoAcid= self.rnaCodonTable[codon] #acess what aa is coded by codon
                self.aaComp[aminoAcid]+=1 #update count of amino acid by 1 if codon is present
                #in this method, if a codon is invalid, it is taken care of by not updating counts
                # Hence, we do not changing the reading frame or delete letters
            
        return self.aaComp

    def cleanedSequence(self,inSeq):
        """
        Cleans up a fasta sequence by unwanted characters and makes all uppercase

        Input: string of nucleotide sequence
        Output: cleaned up, edited sequence

        """
        upperSequence =inSeq.upper() #convert string to upper case

        upperSequence.replace(" ","") #get rid of white space


        emptyPlaceHolder='' #where we will add the edited sequence to

        #Checking if letters are valid nucleotides and adding only valid nucleotides to the emptyPlaceHolder string
        for index in range (0, len(upperSequence)): #for loop to access all indices of string character
            if upperSequence[index] in self.nucComp.keys(): #if character is present as a valid key in self.nucComp Dict
                emptyPlaceHolder+= str(upperSequence[index]) #add the character to the emptyPlaceHolder String

            else:
                emptyPlaceHolder+= "" 

        perfNuc= emptyPlaceHolder #finalized string for us to use
        return perfNuc #this is the final cleaned up sequence where there are no "." "-" and "_" characters
        
        


    
    def nucComposition(self):
        """
        Updates nucleotide composition dictionary with counts of A, T, G ,C, N and U nucleotides
        
        Input: None
        Output: Updated nucComp dictionary with counts
        """

        cleanedSeq=self.cleanedSequence
        for nucleotide in cleanedSeq: 
            if nucleotide in self.nucComp.keys(): #checks if nucleotide is valid. If so, update
                self.nucComp[nucleotide]+=1
        return self.nucComp #return nucleotide composition dictionary
        
    
    def codonComposition(self):
        """
        Updates the codonComp dictionary with counts of valid codons. Uses output of 
        codonListGenerator() method to get a list of codons. 

        input: None
        output: Updated counts of codons in the codonComp Dictionary
        """

        
        codonList= self.codonListGenerator() #generate the list of codons 
    
            
        for codon in codonList:# to update values in codonComp if valid codons exist
            if codon in self.rnaCodonTable.keys():
                self.codonComp[codon]+=1     #only if the codon is a valid codon will we update the dictionary using that specific codon as a key. Invalid codons are not deleted, but do not contribute to the addition in values
        
        return self.codonComp #return codonComp dictionary

    def codonListGenerator(self):
        """
        Returns a list of codons from cleaned Dna Sequence (upper case, no white space)

        input: None
        output: Updated counts of codons

        """
        cleanedSeq=self.cleanedSequence
        rnaSequence='' #define rnaSequence string
        if "T" in cleanedSeq:
            rnaSequence=cleanedSeq.replace("T", "U")  #convert DNA to RNA and add to string by converting T to U
        
        #I will now slice the RNA sequence and create a list containing codons (valid and invalid)
        self.listOfCodons= [rnaSequence[nucleotide:nucleotide+3] for nucleotide in range(0, len(rnaSequence), 3)] #https://stackoverflow.com/questions/43982938/split-string-into-groups-of-3-characters
        
        return self.listOfCodons

    def shannonEntropyCalculator(self,probabilityDictionary):
    	"""
    	Returns a dictionary of the shannon entropy values for each codon position

    	Input: probability dictionary (key(codon number): values(string of codon and its probability))
    	eg. {1:{ATG:1/4, TTT:3/4} 2:{etc}}

    	Output: dictionary of shannon entropy values


    	"""
    	shannonEntropyDict={}
    	for codonNumber, probabiltiesValues in probabilityDictionary.items():
    		for codon, probabilty in probabiltiesValues.items():
    			storageCodon={} #We use this a temporary storage
    			entropy= -probabilty*(math.log2(probabilty)) #calculate the sharnnon entropy for the string of codon
    			storageCodon[codon]=entropy
    			shannonEntropyPosition=sum(storageCodon.values())# Now sum up all the entropy values at the specific codon number
    			shannonEntropyDict[codonNumber]=shannonEntropyPosition #append shannon entropy value for codon number
    		storageCodon.clear()#clear the dictionary for the next codon number :)
    	return shannonEntropyDict

    def probabilityDictCreator(self):
    	"""
    	Creates the probability dictionary for each codon number
		
		probability dictionary (key(codon number): values(string of codon and its probability))
    	eg. {1:{ATG:1/4, TTT:3/4} 2:{etc}}


    	KEY ASSUMPTION:
    	Since the final aligned  sequences from mafft may have the "-" character, we do the following:
    	for a codon position (lets say codon 1) if over 80% of the codons have the "-" character, then we ignore
    	and delete that codon position and its codons and move onto the next set of codons 

    	eg.
    	list of codons before editing:  [[-TAA, -AT, A--], [ATG,TAA,TTT]]

    	after editing we remove the the first index of this list because the fraction of codons in that position that are ambiguous is
    	over 80%

    	Final list: [[ATG,TAA,TTT]]<---- This is now a list containing all the possible codons in position one.

    	"""
    	listsToKeep=[] 
    	prunedCodonListofLists=[]
    	codonPos={}
    	for listCodons in self.codonListofLists:
    		index= self.codonListofLists.index(listCodons)            
    		countOfAmbiguity=0 #keep count of the number of codons that have the alignment "-" character
    		for codon in listCodons:
    			if "-" in codon:
    				countOfAmbiguity+=1


    		fractionOmission=countOfAmbiguity/len(listCodons)
    		if fractionOmission<0.80:
    			listsToKeep.append(index) #we only append indexes if there are fewer than 80% of ambigious codons for codons at that codon position

    	for index in listsToKeep:
    		prunedCodonListofLists.append(self.codonListofLists[index]) #append into a list of lists

    	for codonList in prunedCodonListofLists:
    		position=(prunedCodonListofLists.index(codonList))+1 #the position of the codon is the index +1
    		codonPos[position]=codonList

    	# now we calculate the probabilites of all codons at each codon position
    	probabilityDictionary={key:{} for key in codonPos}
    	for codonPosition,listOfCodons in codonPos.items():
    		for codon in listOfCodons:
    			codonFreq=listOfCodons.count(codon)
    			probability=codonFreq/len(listOfCodons)
    			probabilityDictionary[codonPosition][codon]=probability
    	return probabilityDictionary







       
    
    def nucCount(self):
        """ 
        Returns the total nucleotide count by summing values of the nucComp Dictionary
        input: none 
        output: count of the sum of all nucleotide frquencies
        """
        
        self.nucComposition() #call method for nucComp that will return updates sequence
        
        #now I will add the sum of all the nucleotide values in the nucComp Dictionary
        self.nucCountValue=sum(self.nucComp.values()) #https://stackoverflow.com/questions/4880960/how-to-sum-all-the-values-in-a-dictionary
        return self.nucCountValue

    def codonTypesList(self,inSeq):
        """
        Determines the nth codon of a string of dna 
        
        input: string of DNA sequence
        output: dictionaries for codons and amino acids. key is an int representing the nth codon. Value is the string of aa or codon
        """
        codons = [inSeq[i:i+3] for i in range(0, len(inSeq), 3)] #codon is every 3 nucleotides

        #now we determine the codon for each position in CDS
        
        if self.counter==0:
            for codon in codons:
                self.codonListofLists.append([codon])
                
            self.counter+=1
        else:
            codonNumber=0
            for codon in codons:    
                self.codonListofLists[codonNumber].append(codon)
                
                codonNumber+=1

        

        
                




                

        return self.codonListofLists


class ProteinParam:
    """
    This class contains methods to compute the pI, molecular weight, composition and 
    molar extinction, mass extinction of an amino acid sequence
    
    Input: String of amino acid sequence
    Output: Object that can compute pI, molecular weight, composition and 
    molar extinction, mass extinction of an amino acid sequence
    
    class attributes:
    molecular weight dictionary 
    pKa dictionary for positively charged amino acids
    pkA dictionary for negatively charged amino acids
    aaNterm pKa
    aaCterm pkA
    
    """
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        """ 
        initializes object and creates aaComposition dictionary as an attribute of class
        parses the input string immediately to make upper case and calculate counts of each
        amino acid in the dictionary
        
        input: string of the protein sequence
        output: object of class the updated dictionary for amino acid composition
        
        """
        
        #creates a dictionary of the amino acids from aa2w and sets count as 0. Updates values in dict +1 for each occurence of letter 
        self.protein = protein #initlizes string as attribute of class
        parsedProtein = protein.replace(" ","") #edits sequence to remove space
        upperProtein= parsedProtein.upper() #casts characters into upper case
        self.aaCompositionDict= {key:0 for key,value in self.aa2mw.items()} #citation below, but this creates new dictionary using keys from aa2mw dictionary
        #sets value (count of each aa) to 0 first
        #citation:https://stackoverflow.com/questions/12117080/how-do-i-create-dictionary-from-another-dictionary
        for aa in upperProtein: #for loop to scan through the cleaned protein sequence
            if aa in self.aaCompositionDict.keys(): 
                  self.aaCompositionDict[aa]+=1  #updates dictionary with counts if amino aa is a key in the aacomp dict
        

    def aaCount (self):
        """
        Computes the counts of the aminio acids by summing up value from the aaComposition Dictionary
        
        Input: none
        output: aaCounts is the sum of all amino acid counts
        """
        self.aaCounts=0 #sets counting variable (an object of class int)
        for value in self.aaCompositionDict.values(): #loop through the values of the aaComposition Dict
            self.aaCounts+= int(value) #add to the counting variable
        return self.aaCounts #returns the counts

    def pI (self):
        """
        Calculates pI, the pH that yields lowest charge. This method uses a search for best method
        where we iterate through all pH values from 0.00 to 14.00 until we hit the pH that 
        produces the lowerst charge. Returns 0 if self.aaCounts is 0, meaning no valid amino acids in
        user input. 

        Input: none
        Output: pI value (pH that produces lowest charge)
        
        """
        
       # I got the following code from April 17th lecture from BME 160
        bestCharge= 9999999999999999 #set best charge as high as possible

        for pH100 in range(0,1400+1): #iterate over all possible numbers from 0-1400
            pH= pH100 * 10**(-2) #we thus can get pH to 2 decimals in range 0.00-14.00 because of above
            thisCharge= (self._charge_(pH)) #generate charge using pH as parameter 
            if abs(thisCharge) < bestCharge: #calculate abs value of charge as charge can be negative
                bestCharge= thisCharge 
                desiredpH= pH #this is the pI
        if self.aaCounts>0:
            return desiredpH
        else: #if no valid AA in input string, return a 0
            return 0


        
    
    def aaComposition (self) :
        """
        Returns the aaCompositionDict generated from initialize
        input: none
        Output: Returns aaCompositionDict
        
        """
        return self.aaCompositionDict

    def _charge_ (self,pH):
        """
        Calculates charge of protein sequence if given a parameter of pH
        
        input: pH
        output: net charge 
        
        """
        negativeAACharge=0 #sets initial charge of negative amino acids as 0
        positiveAACharge=0 #sets initial charge of positively charged amino acids as 0
        
        #Logic of method is for the amino acid of the key in the Neg Dict, we access the
        #amino acid as a key in the aacomp dictionary and the aa2chargeNeg dictionary
        for aa in self.aa2chargeNeg.keys():
            negativeAACharge+= (self.aaCompositionDict[aa]*(10**pH/((10**self.aa2chargeNeg[aa])+(10**pH)))) #returns charge of negative amino acids
        
        for aa in self.aa2chargePos.keys():
            positiveAACharge += (self.aaCompositionDict[aa]*(10**self.aa2chargePos[aa]/((10**self.aa2chargePos[aa])+(10**pH)))) #returns charge of positive amino acids
            
        NtermCharge= (10**self.aaNterm/(((10**self.aaNterm) + (10**pH)))) #calculates Nterminus charge
        CtermCharge= (10**pH/(((10**self.aaCterm) + (10**pH)))) #calculates the Cterminus 
                               
        
        negativeCharge= negativeAACharge + CtermCharge #calculates the total negative charge
        positiveCharge= positiveAACharge + NtermCharge #calculates the total positive charge
        
        totalCharge= positiveCharge-negativeCharge #total charge is the difference between postive charge- negative charge
        return totalCharge

    def molarExtinction (self,cystine=True):
        """
        Calculates molar extinction under a default parameter of oxidizing conditions of cystine
        
        Input (default parameter): cystine= True (oxidizing conditions
        Output: molar extinction value
        """
        if cystine==True: #calculates molar extinction under oxidizing conditions
            cystineExtinction = self.aa2abs280['C']*self.aaCompositionDict['C'] #access values using same key (ie. multiple count of aa by the extinction value)
            tyrosineExtinction= self.aa2abs280['Y']*self.aaCompositionDict['Y']
            trptophanExtinction= self.aa2abs280['W']*self.aaCompositionDict['W']
            sum= cystineExtinction + tyrosineExtinction + trptophanExtinction
        else: #when cystine is false (ie. reducing conditions), cysteine does not contribue to extinction
            tyrosineExtinction= self.aa2abs280['Y']*self.aaCompositionDict['Y']
            trptophanExtinction= self.aa2abs280['W']*self.aaCompositionDict['W']
            sum= tyrosineExtinction + trptophanExtinction
        return sum
    
    

    def massExtinction (self,cystine=True):
        """
        Calculates the mass extinction coefficient
        input (default parameter): cystine= True (oxidizing conditons)
        Ouput: mass extinction value which is the molar extinction/molec weight of peptide
        """
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        """
        Function calculates number of peptide bonds and water produced and subtracts it from weight of aa to get
        mol weight of peptide
        
        input: none
        output: molecular weight of peptide
        """
         
        
        numberOfPeptideBonds= int(self.aaCount()-1) # calculate number of peptide bonds. This is equivalent to moles of water
        molWeightWaterInPeptide= self.mwH2O*numberOfPeptideBonds #calculate the molecular weight of water from peptide
        molWeightofAminoAcids= sum(self.aaCompositionDict[k]*self.aa2mw[k] for k in self.aa2mw) #https://stackoverflow.com/questions/16087118/multiplying-and-then-summing-values-from-two-dictionaries-prices-stock
        molecularWeightPeptide= molWeightofAminoAcids-molWeightWaterInPeptide
        if self.aaCount()>0:  #returns molecular weight of peptide if we have at least 1 valid aa
            return molecularWeightPeptide
        else:  #handles case if no valid amino acids are given in string input.
            return 0


class FastAreader :
	''' 
	Define objects to read FastA files.
	
	instantiation: 
	thisReader = FastAreader ('testTiny.fa')
	usage:
	for head, seq in thisReader.readFasta():
		print (head,seq)
	'''
	def __init__ (self, fname=''):
		'''contructor: saves attribute fname '''
		self.fname = fname
			
	def doOpen (self):
		''' Handle file opens, allowing STDIN.'''
		if self.fname is '':
			return sys.stdin
		else:
			return open(self.fname)
		
	def readFasta (self):
		''' Read an entire FastA record and return the sequence header/sequence'''
		header = ''
		sequence = ''
		
		with self.doOpen() as fileH:
			
			header = ''
			sequence = ''
			# skip to first fasta header
			line = fileH.readline()
			while not line.startswith('>') :
				line = fileH.readline()
			header = line[1:].rstrip()

			for line in fileH:
				if line.startswith ('>'):
					yield header,sequence
					header = line[1:].rstrip()
					sequence = ''
				else :
					sequence += ''.join(line.rstrip().split()).upper()

		yield header,sequence













