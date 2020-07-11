
#!/usr/bin/env python3
# Name: Nitin Vidyasagar (nvidyasa)
# Group Members: Hannah Knotter
"""
Purpose:
This program will take in a normalized pickled data frame of gene counts and plot a scatter plot 
in log scale for different libraries

Required modules:
pandas as pd
seaborn as sns
matplotlib.pyplot as plt
csv 

run as:

python3 geneCountsPlot.py inFile1.pkl inFile2.csv inFile3.csv inFile4.csv outPrefix

"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
import csv


class CommandLine():
    """Handle the command line to determine inFiles and outFile.
    Only arguments that are specified in the command line are the file names."""
    def __init__(self, inOpts=None):
        """Implement a parser to interpret command line argv string using argparse."""
        import argparse
        self.parser = argparse.ArgumentParser(description = 'parses csv and *.jam files from assignReadsToGenes5.py', 
                                             
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = 'python3 %(prog)s <input >output'
                                             )
        self.parser = argparse.ArgumentParser(description="parses pickled files of data frames and generates JPEG scatterplot")
        self.parser.add_argument("inFile1", action="store", nargs="+", help="pickled gene data frame")
        self.parser.add_argument("inFile2", action="store", nargs="+", help="csv file for Genes we want to highlight")
        self.parser.add_argument("inFile3", action="store", nargs="+", help="csv file for colors of Genes")
        self.parser.add_argument("inFile4", action="store", nargs="+", help="csv file library genotypes")
        self.parser.add_argument('-x', '--xAxis', action = 'store', nargs='?',help='determines what genotype from library to plot on x axis')
        self.parser.add_argument('-y', '--yAxis', action = 'store', nargs='?',help='determines what genotype from library to plot on y axis')
        self.parser.add_argument("outPrefix", action="store", nargs="+", help="outPrefix For Pdf file")
        self.parser.add_argument('-a', '--addOne', action = 'store', nargs='?', const=True, default=False, help='adds one to all values in data if selected to handle counts of 0')

    

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class CSVtoDictionary:
    """ 
    This class will create an object that can parse a csv file and return a dictionary 

    to initialize:
    object=  CSVtoDictionary()

    Methods include:
    1) csvSingleDictionary This method creates a dictionary for that has the key as the first item in each row of csv. Value is second
    item in the row

    2)csvListDictionary:  This method creates a dictionary for that has the key as the header  of a csv. Values is items in the same column

    """


       
    def csvSingleDictionary(self,inFile): 
        """
        Parses a csvFile for where key is the the first time of a row and value is the second item of the row

        input: inFile(this should be csv file)
        Output: Dictionary for first value in row is key and value is the second value of row

        RELEVANCE:
        This is to be used to generate the genotype translations of the sequencing libraries (for commandline) and also 
        specifies what color we are to use for each gene we want to highlight 

        """
        #citation:
        #https://stackoverflow.com/questions/46943126/creating-a-dictionary-from-csv-using-python


        f= open(inFile, 'r')   #open csv file

        reader= csv.reader(f) #stor output of reader method from csv module

        dictStore= {} #this is where we will append the key value pairs

        for row in reader: 
            dictStore[row[0]]=row[1].replace(" ","") #clean up any white spaces from entries

        return dictStore

       
    def csvListDictionary(self,inFile):   
        """
        Parses a csvFile for where key is the column name and values is a list of entries in the column

        input: inFile(this should be csv file)
        Output: Dictionary for key is the column name and values is a list of entries in the column

        RELEVANCE:
        This is used to mark what our genes of interest are that we want to highlight on a plot. 

        """

        #citation:https://stackoverflow.com/questions/56141590/how-to-read-csv-file-into-dataframe-using-pandas

      
        data = pd.read_csv(inFile)
     
        
        dictionaryStorage = {col: list(data[col]) for col in data.columns}

        for key, value in dictionaryStorage.items():
            dictionaryStorage[key]= [x.replace(" ","") for x in value if str(x) != 'nan']  #https://stackoverflow.com/questions/57888161/how-can-i-remove-nan-values-from-some-columns-in-a-csv-file
           

      
        return dictionaryStorage


class PickleToDataFrame:
#     """ Creates an object that can open a pickled data frame and access contents """
	def __init__(self,inFile1):
#         """Initializes inFiles for access in class methods """
		self.inFile1= inFile1


	def dataFrameGenerator(self):
		"""Reads a pickled dataframe to return a pandas dataframe"""
		self.dataFrame= pd.read_pickle(self.inFile1)
		return self.dataFrame


class Scatterplot:
    """
    Creates an object that can plot a log normalized scatterplot using seaborn
    
    to initalize:
    object= Scatterplot(df,favGeneDict,colorDict,libraryDict,x,y,outPrefix)

    This class contains methods:

    1) dataFrameLibraryConversion: This method will rename columns of a dataframe using inputs from a csv file to specify genotypes

    2) metaGeneSum: This calculates the sum of all proteins of the same gene class (specified in csv) and thus collapses normalized readcounts
    of genes in the same group to create a metagene. it also drops and removes from the data frame values that we collapsed into
    the metagene.

    3)plotGenerator: This method generates the scatterplot of the two libraries/genotypes we specified on the commandline. genes of
    interest (marked using csv) are highlighted on the plot.


    """

    def __init__(self,df,favGeneDict,colorDict,libraryDict,x,y,outPrefix):
        """Initializes the class attributes so we can access them in class methods"""
        self.df = df #dataframe containing RPM normalized genecounts
        self.favGeneDict= favGeneDict #this contains the genes arranged according to groups. Key is geneGroup and value is list of genes
        self.colorDict= colorDict #This specifies what color to label the geneGroup. Key is geneGroup, value is colors in seaborn
        self.libraryDict= libraryDict #This contains what is the genotype of the libraries. 
        self.x=x #what genotype will be plotted as x axis
        self.y=y #what genotype will be plotted as y axis
        self.outPrefix=outPrefix #this is output file name for us to use. 

    def dataFrameLibraryConversion(self):
        """
        This method will change the library names to the necessay genotypes

        Input: none
        Output: DataFrame with columns renamed by mapping names to the library dict 

        """
        self.new= self.df.rename(columns = self.libraryDict) #we re-named columns byt mapping the data frame columns to the values in the dictionary we generated from the csv
    
        return self.new

    def metaGeneSum(self):
        """
        This calculates the sum of all proteins of the same gene class (specified in csv) and thus collapses normalized readcounts
        of genes in the same group to create a metagene. it also drops and removes from the data frame values that we collapsed into
        the metagene.

        The approach is to add a column called Protein where we add the key from the favGeneDict to all the genes in the values. Then
        we utilitze the groupby and sum appoach to sum all the values.


        Input: none
        Output: self.edited which the new dataframe with the metagene appended on the original dataframe. Genes present in faveGeneDict are dropped

        """

        dfOld= self.new.copy() #we make a copy of this so as to append the metagene data frame onto it

        for wbGeneID, normCount in self.new.iterrows(): #iterate through rows of data frame
            for geneClass, listGenes in self.favGeneDict.items(): #iterate through faveGenedict
                if wbGeneID in self.favGeneDict[geneClass]: 
                    self.new.loc[wbGeneID, 'Protein']= str(geneClass) #label gene with the gene class it belongs to in the newly created protein column

        self.new=self.new.groupby(['Protein']).sum() #sum up all the genes in the metagene
        self.new['Protein'] = self.new.index #make the metagene names as an index so we can append to the original data frame
        self.new['Size']='Big' #For those in the metagene/fav gene dataframe, I mark the size as big so we can increase the size of points when we plot

        pd.set_option("display.max_rows", None, "display.max_columns", None)
        print(self.new)
        

        self.new=self.new.append(dfOld,sort=False) #we append the old, original dataframe to the new data frame
       
       
        self.edited= self.new.fillna('other') #we will all NaN values with other at this point


    
        return self.edited

    def plotGenerator(self,addArg):
        """
        Generates a scatter plot using seaborn. The x and y axis are parsed in from command line as class attributes. The output will be saved
        as a jpg image of the graph

        Since we log transform, there are 2 approaches to the data: 
        1) We add 1 to all values in the data to remove any 0's   (IF ADDARG IS TRUE)
        2) We remove indexes from data frame that have value of 0 (IF ADDARG IS FALSE)

        The choice of either is passed in from commandline

        Input: addArg (boolean)
        Output: jpg of graph


        """

        if addArg==True:

            #https://stackoverflow.com/questions/30794525/adding-one-to-all-the-values-in-a-dataframe

            numeric_cols = [col for col in self.edited if self.edited[col].dtype.kind != 'O'] #Basically gets the list of values in df that are not objects
            self.edited[numeric_cols] += 1 #add one to every value in the data frame


      


            plt.figure(figsize=(5,5)) #we set the size of the figure

            #we will plot a scatterplot and we get rid of the legend, set the color of dots according to color dict and specifying the size of the points
            ax=sns.scatterplot(x=self.x, y=self.y, data=self.edited,
                        hue='Protein', size='Size', sizes= (0.5,40),palette= self.colorDict, legend=False,linewidth=0,edgecolor=None)


            ax.set(xscale="log", yscale="log")
            

            #I am annotating the dots I showed by color. The approach is to keep the text at the same x position and vary the y positon by the same amounts ie 
            #decrease by half.
            yPos=40000
            for key, value in self.colorDict.items():
                if key== "other":
                    pass
                else:
                    ax.text(10**-0.9,yPos, key, color=self.colorDict[key],fontsize=6)
                    yPos-=yPos/2


            
            

            plt.title('15-18nt Counts (RPM) Normalized', fontsize=12)

            #we set the limits for the graphs 
            plt.ylim(10**-1, 10**5)
            plt.xlim(10**-1, 10**5)
            plt.ylabel(u'Δ '+self.y) #I add delta to show delete symbol 

            plt.savefig(self.outPrefix + ".jpg") 


        if addArg==False:

            
            #https://stackoverflow.com/questions/49841989/python-drop-value-0-row-in-specific-columns
            #in this code we drop any values that are 0 because the log plot cannot handle that
            self.edited=self.edited.loc[(self.edited[[self.x, self.y]] != 0).all(axis=1)]

      
            plt.figure(figsize=(5,5))#we set the size of the figure
            ax=sns.scatterplot(x=self.x, y=self.y, data=self.edited,
                        hue='Protein', size='Size', sizes= (0.5,40),palette= self.colorDict, legend=False,linewidth=0,edgecolor=None)


            ax.set(xscale="log", yscale="log") #set log scale
            plt.ylabel(u'Δ '+self.y) #I add delta to show delete symbol 
            plt.xlabel(u'Δ '+self.x) #I add delta to show delete symbol 


            yPos=40000
            for key, value in self.colorDict.items():
                if key== "other":
                    pass
                else:
                    ax.text(10**-0.9,yPos, key, color=self.colorDict[key],fontsize=6)
                    yPos-=yPos/2


            plt.title('15-18nt Counts (RPM) Normalized', fontsize=12)

           
            plt.ylim(10**-1, 10**5)
            plt.xlim(10**-1, 10**5)


            plt.savefig(self.outPrefix + ".JPEG") 

    def logFoldChangeCalculator(self):
        """
        This method calculates the log2 fold change for the genes of interest. To do this, we iterate through our dictionary and
        find the value at the edited data frame  for the x axis and y axis library. We then normalize using a log2 
        function to figure out the fold change. The information on the fold change will be stored as a list of lists

        input: none
        Output: List containing the geneclass and the foldchange compared to the x axis

        """
        foldChangeList=[] #this list will store information on the geneClass and the log2 fold change we calculate
        
        for geneClass, listGenes in self.favGeneDict.items():
            xRPM= self.edited.loc[str(geneClass), str(self.x)] #access spcific RPM value using index and column (citation# https://pythonhow.com/accessing-dataframe-columns-rows-and-cells/)
            yRPM= self.edited.loc[geneClass,self.y] #access spcific RPM value using index and column
            foldChange=math.log2(yRPM/xRPM) 
            foldChangeList.append([geneClass,foldChange])
        return foldChangeList



def main(inCL=None):
    """
    This function parses information from commandline and files. The function also prints out into the screen the 
    log2 fold change values the genes of interest

    Input: None
    Output: Printed statements on screen for fold change

    """
    
    if inCL is None:
        myCmdLine = CommandLine()
    else:
        myCmdLine = CommandLine(inCL)

    #unpickle my pickle objects. Pass in dataframe from commandline as argument
    geneCountsPickle=PickleToDataFrame(myCmdLine.args.inFile1[0])
    

    #creates object of CSVtoDictionary Class
    csvParse= CSVtoDictionary() 

    #create our data frame from pickled dataframe
    dataFrame=geneCountsPickle.dataFrameGenerator()
    #create variables from commandline arguments
    y=myCmdLine.args.yAxis
    x=myCmdLine.args.xAxis
    outPrefix=myCmdLine.args.outPrefix[0]
    addArg= bool(myCmdLine.args.addOne)
    

    #create the dictionaries from CSV files
    favGeneDict= csvParse.csvListDictionary(myCmdLine.args.inFile2[0]) #this will be a dictionary specifying the genes we want to highlight and collapse into a metagene
    colorDict=csvParse.csvSingleDictionary(myCmdLine.args.inFile3[0]) #this will be a dictionary specifying the color for the genes of interest
    libraryDict= csvParse.csvSingleDictionary(myCmdLine.args.inFile4[0]) #this will be a dictionary specifying the library and the genotype associated with it

  
    #Create object of Scatterplot Class
    sca= Scatterplot(dataFrame,favGeneDict,colorDict,libraryDict,x,y,outPrefix)
   
    new= sca.dataFrameLibraryConversion()


    metaGeneDf= sca.metaGeneSum()
    print(metaGeneDf.head())
    

    foldList=sca.logFoldChangeCalculator() #create the fold list for the genes of interest
    for item in foldList:
        geneClass=item[0]
        fold=item[1]
        print("The log2 foldchange between for " + str(geneClass) +" is {:.2f}".format(fold))
    
    sca.plotGenerator(addArg)

main()










