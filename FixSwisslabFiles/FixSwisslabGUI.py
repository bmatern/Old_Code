import sys
import os
#import numpy
#import subprocess
import getopt
import random
#from subprocess import check_output
from os.path import split, join


import Tkinter, Tkconstants, tkFileDialog
from Tkinter import *


#import random
#import math
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC

#from gc import garbage
#from Bio import SeqRecord
#import operator

class FixSwisslabGUI(Tkinter.Frame):

    # Initialize the GUI
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("Fix the Swisslab Export Files")
        self.parent = root
        self.initialize(root)
        
    # Initialize GUI elements
    def initialize(self, root):

        # This is the directory the python executable is running from.
        FileAndPath = os.path.abspath(__file__)
        self.idir, self.ifile = os.path.split(FileAndPath)

        # GUI options
        self.label_opt = {'fill': Tkconstants.BOTH, 'padx': 10, 'pady': 6}
        self.button_opt = {'fill': Tkconstants.BOTH, 'padx': 50, 'pady': 6}
        self.dir_opt = {'initialdir': self.idir,  
            'parent': root, 'title': 'Choose an input file....'}

        self.master = root
        self.frame = Tkinter.Frame(self.master)

        self.instructionText = Tkinter.StringVar()
        self.instructionText.set('Select an input file.')
        #self.defaultMessage()
        Tkinter.Label(self.frame, width=100, height=1, textvariable=self.instructionText).pack(self.label_opt)
        
        self.submitButton=Tkinter.Button(self.frame, text='Choose Input File', command=self.getInputFile)
        self.submitButton.focus_set()
        self.submitButton.pack(self.button_opt)
        
        self.inputFileText = Tkinter.StringVar()
        self.inputFileText.set('Broken Swisslab file:')        
        Tkinter.Label(self.frame, width=100, height=1, textvariable=self.inputFileText).pack(self.label_opt)
        
        self.inputDirectoryStringvar = Tkinter.StringVar()
        Tkinter.Entry(self.frame, width=100, textvariable=self.inputDirectoryStringvar).pack(self.label_opt)
  
        self.outputFileText = Tkinter.StringVar()
        self.outputFileText.set('I will put the fixed file here:')        
        Tkinter.Label(self.frame, width=100, height=1, textvariable=self.outputFileText).pack(self.label_opt)
        
        self.outputDirectoryStringvar = Tkinter.StringVar()
        Tkinter.Entry(self.frame, width=100, textvariable=self.outputDirectoryStringvar).pack(self.label_opt)
                
        self.fixButton=Tkinter.Button(self.frame, text='Fix the file', command=self.fixFileClick)
        #self.submitButton.focus_set()
        self.fixButton.pack(self.button_opt)
        

        # This is what happens when user hits "Enter" on the keyboard.  Simple stuff.
        def returnEvent(event):
            print('Returned')
            self.getInputFile()
            
        self.master.bind('<Return>', returnEvent)

        self.submitButton.pack(self.button_opt)

        self.frame.pack()
        
    # getInputFile method is called when the user presses the input directory button
    def getInputFile(self):
        self.submitButton.config(state='disabled')
        self.setInputDir(tkFileDialog.askopenfilename(**self.dir_opt))
        #self.setInputDir(tkFileDialog.askdirectory(**self.dir_opt))
        self.submitButton.config(state='normal')
        
    def fixFileClick(self):
        inputFileName = self.inputDirectoryStringvar.get()
        outputFileName = self.outputDirectoryStringvar.get()
        
        print ('Im gonna fix these files:\n' + inputFileName + '\n' + outputFileName)
        
        fixFile(inputFileName, outputFileName)
        
        

    # setInputDir is a subprocess to be called when an input directory is selected.
    # Look into the directory contents, and present an extractor window to the user.
    def setInputDir(self, inputFile):
        inputFile = os.path.normpath(inputFile)
        
        self.inputDirectoryStringvar.set(inputFile)
        
        inputFileDir, inputFileName = os.path.split(inputFile)
        
        outputFileDir = os.path.join(inputFileDir, 'FixedFiles')
        outputFileName = inputFileName.replace('.csv','_fixed.csv')
        
        self.outputDirectoryStringvar.set(os.path.join(
            outputFileDir, outputFileName))
        #Set the text variables in here.

    
        
        
        
def fixFile(inputFileName, outputFileName):
    print ('Time to fix this file: ' + inputFileName)
    
    #Read the input.  I don't need to read and process in separate steps, but might as well.
    inputTextLines = []
    outputTextLines = []
    
    try:
        fileInputObject = open(inputFileName, 'r')
        for line in fileInputObject:
            inputTextLines.append(line)        
        fileInputObject.close()
    except Exception:
        print('I couldn\'t open the file.  Does it exist?  Permission issue?')
        raise
    
    if (len(inputTextLines) < 1): 
        print ('This file looks empty.  That\'s not normal.')
        raise Exception('Empty Input File')
    
    for index, line in enumerate(inputTextLines):
        
        
        print('Line Index:' + str(index) + '\n')
        print('Before Change:\n' + line + '\n')
        
        newLine = ''

        #print(line)
        
        # Skip processing of the first two lines.  
        # Maybe I'll need to do something with the header in here.
        if(index == 0):
            # Add the Locus header
            newLine = line.replace(',\"NC1\"',',\"Locus\",\"NC1\"')
        elif(index == 1):
            newLine = line
        else:
            # They put a comma in the middle of some of the data columns.  
            # So I need to split by ","
            stringTokens = line.split('","')
            #newLineTokens = []
            
            #print('This line has ' + str(len(stringTokens)) + ' columns.')
            
            # The locus information is in the 14th column (13 is zero based)
            # It is before the area where :=: appears.
            alleleColumn = stringTokens[13].split(':=:')[0]
            # There are commas in this (It's a csv file :( , i don't know what they were thinking)
            # I split by space later, so i'll put a space there instead.
            alleleColumn = alleleColumn.strip().replace(',',' ')
            #print ('The allele information:' + alleleColumn)
            
            currentGeneLoci = []
            
            # Each allele is separated by a space
            for alleleText in alleleColumn.split(' '):
                #print ('Found this allele:' + alleleText) 
                
                alleleSplit = alleleText.split('*')
                
                if(len(alleleSplit) == 2):
                    currentGeneLocus = alleleSplit[0]
                    #print ('Found this locus:' + currentGeneLocus)
                    currentGeneLoci.append(currentGeneLocus)
                    
                    #
                else:
                    pass
                    #print('Something wrong with this allele:' + alleleText)     
                    
            #Eliminate Duplicate Genes by converting it to a set and back
            currentGeneLoci = sorted(list(set(currentGeneLoci)))
            
            print('I found a total of ' + str(len(currentGeneLoci)) + ' loci.')
            
            for columnIndex in range (0,3):
                newLine = newLine + (stringTokens[columnIndex] + '","')
                
            # Add the genes we just found to the new data file.
            # If there aren't any, we need to leave this column blank.
            if(len(currentGeneLoci) == 0):
                newLine = newLine + '","'
            else:
                for newLocus in currentGeneLoci:
                    print('appending locus:' + newLocus)
                    newLine = newLine + (newLocus) + ','            
                newLine = newLine[0:len(newLine) -1] + '","'
            
            for columnIndex in range(3,len(stringTokens)):
                newLine = newLine + (stringTokens[columnIndex])
                if(columnIndex < len(stringTokens) - 1):
                    newLine = newLine + '","'
            
            #for dataColumn in stringTokens:
            #   print ('\t' + dataColumn)
            
        print('Edited Line:\n' + newLine + '\n')
        
        outputTextLines.append(newLine)
        
    # Here I am done editing the lines.  
    # Time to write to an output file.
    
    try:
        outputFileObject = createOutputFile(outputFileName)
        
        for outputLine  in outputTextLines:
            outputFileObject.write(outputLine)    
        outputFileObject.close()
        
    except Exception:
        print('I couldn\'t write to the output file.  Does it exist?  Permission issue?')
        raise
    
# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
        tempDir, tempFilename = split(outputfileName)
        if not os.path.isdir(tempDir):
            os.mkdir(tempDir)
        resultsOutput = open(outputfileName, 'w')
        return resultsOutput


