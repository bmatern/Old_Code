import sys
#import os
#import numpy
#import subprocess
import getopt
#import random
#from subprocess import check_output
#from os.path import split, join

import Tkinter

from FixSwisslabGUI import FixSwisslabGUI
#import random
#import math
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC

#from gc import garbage
#from Bio import SeqRecord
#import operator

SoftwareVersion = "FixSwisslabFiles Version 1.0"

def usage():
    print("Oops, you are using this program incorrectly.  See Ben or Mathijs for help.")    



# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
    global inputFileName
    global outputFileName


    inputFileName  = None
    outputFileName = None


    if(len(sys.argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        usage()
        #raise Exception('Problem with Commandline Arguments')
        return False    

    # getopt.getopt(..) is a function for parsing args the way smart people do it
    # For More info, use google or 
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"hvi:o:"
            ,["help", "version", "inputFile=","outputFile="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            elif opt in ("-i", "--inputFile"):
                inputFileName = arg
            elif opt in ("-o", "--outputFile"):
                outputFileName = arg
                
            else:
                print('Unknown commandline option: ' + str(opt) + ':' + str(arg))
                #raise Exception('Problem with Commandline Arguments')

    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    print('Input File:' + str(inputFileName))
    print('Output File:' + str(outputFileName))

    # Quick sanity check.
    if(len(inputFileName) < 4):
        print('Input file is too short:' + str(inputFileName))
        return False
    if(len(outputFileName) < 4):
        print('Output file is too short:' + str(outputFileName))
        return False
    if(inputFileName == outputFileName):
        print('The input and output file are the same.  I\'m not comfortable with that, please change the output to something new.')
        return False


    return True

    

if __name__=='__main__':
    try:

        if(not readArgs()):
            print('Something wrong with the commandline arguments.  Let\'s just launch the GUI.')
            
            
            root = Tkinter.Tk()
            FixSwisslabGUI(root).pack()
            root.mainloop()
            #raise Exception('Problem with Commandline Arguments')
        else:
            FixSwisslabGUI.fixFile(inputFileName, outputFileName)
        
        
        
        print('Done.  Ben is an excellent programmer.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise


