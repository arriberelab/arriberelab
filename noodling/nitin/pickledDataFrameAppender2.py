#!/usr/bin/env python3
# Name: Nitin Vidyasagar (nvidyasa)
# Group Members: Hannah Knotter



"""
200426

Purpose: 
To parse pickled dataFrames of same dimensions and append them together to create a 
master data frame of the same dimensions. Master data frame is then pickled for future analysis

run as:

python3 pickledDataFrameAppender.py outFile.pkl 

Required modules:
import sys, csv, time
import pandas as pd
import pickle
from logJosh import Tee
import os

"""
import sys, csv, time
import pandas as pd
import pickle
import os

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
		self.parser = argparse.ArgumentParser(description="concatenates dataframes along axis=1")
		self.parser.add_argument("-p",'--directoryPath', action = 'store', nargs='?', help='provide string path for dataframes to append in same directory. Make sure all dataframes are in one directory and only dataframes are in that directory')
		self.parser.add_argument("outFile", action="store", nargs="+", help="output file is pickled data frame: outFile.pkl")


		if inOpts is None:
			self.args = self.parser.parse_args()
		else:
			self.args = self.parser.parse_args(inOpts)

class PickleToDataFrame:
	""" 
	Creates an object that can open a pickled data frame and access content
	
	To instantiate:

	object= PickleToDataFrame()
	"""

	def dataFrameGenerator(self,path):
		"""Reads a pickled dataframes from path to return a master pandas dataframe but concatenating the values"""
		#https://stackoverflow.com/questions/32444138/concatenate-a-list-of-pandas-dataframes-together
		listOfDfs=[] #store data frames as a list
		files = os.listdir(path)
		for file in files: #we loop through dataframes in that directory and then create a 
		    df=pd.read_pickle(file)
		    listOfDfs.append(df)

		eDataFrame=pd.concat(listOfDfs, axis=1)

		#Convert ALL NaN Values as 0
		self.editedData=eDataFrame.fillna(0.0)
		
		return self.editedData


def main(inCL=None):
	"""
	Appends data frames together and saves the master dataframe as pickled dataframe.

	input: none
	output: pickled dataframe containing concatenated counts of all libraries


	"""
	if inCL is None:
		myCmdLine = CommandLine()
	else:
		myCmdLine = CommandLine(inCL)

	directoryPath=myCmdLine.args.directoryPath
	print(directoryPath)
	

	pickledDataFrameAppender= PickleToDataFrame()

	masterDataFrame= pickledDataFrameAppender.dataFrameGenerator(directoryPath)

	print(masterDataFrame)

	return masterDataFrame.to_pickle(myCmdLine.args.outFile[0]) 


main()
