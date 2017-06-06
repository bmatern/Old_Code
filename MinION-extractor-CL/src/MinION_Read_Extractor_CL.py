# This file is part of MinION-extractor-CL.
#
# MinION-extractor-CL is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MinION-extractor-CL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MinION-extractor-CL. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0

SoftwareVersion = "MinION-extractor-CL Version 1.0"

import getopt
import os
from os import listdir
from os.path import isfile, join, split
import sys
from MinION_Read_Extractor import *



# This script can be used to extract the fastq and fasta sequences 
# from a directory of fast5 files, from MinION Output
# It expects a single parameter, to specify a search directory.  
# Output files are generated in a read_extracts folder under the input directory.


# Read Barcode SampleNames
def readBarcodeSampleNames(barcodeFileNameWithPath):
    try:
        #runningDirectory = os.path.dirname(sys.argv[0])
        #barcodeFileNameWithPath = join(runningDirectory, 'barcodes.txt')
        # print ('running directory:' + runningDirectory)

        barcodes = {}
        barcodeFile = open(barcodeFileNameWithPath, 'r')
        for i, line in enumerate(barcodeFile):
            if (len(line.strip()) > 0):
                # Filter comments.
                if not (line.strip()[0:1]=='#'):
                    tokens = line.split()
                    barcodes[tokens[0]] = tokens[1]

        barcodeFile.close()
        return barcodes

    except Exception, e:
        print ('Problem reading barcode sample file: ' + str(e))
        raise


# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
    global inputRootDirectory
    global outputResultDirectory
    global minimumReadLength
    global maximumReadLength
    global minimumQuality
    global barcodeSampleMapFilename
    global barcodeSampleMap
    global rundate
    
    inputRootDirectory       = None
    outputResultDirectory    = None
    minimumReadLength        = 0
    maximumReadLength        = 0
    minimumQuality           = 0
    barcodeSampleMapFilename = None
    barcodeSampleMap         = None
    rundate                  = None
    

    if(len(sys.argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        usage()
        return False    


    # getopt.getopt(..) is a function for parsing args the way smart people do it
    # For More info, use google or 
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"m:M:hvi:o:b:r:q:"
            ,["minlen=", "maxlen=", "help", "version", "idir=","odir=","barcodes=", "rundate=", "minqual="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            elif opt in ("-i", "--idir"):
                inputRootDirectory = arg
            elif opt in ("-o", "--odir"):
                outputResultDirectory = arg
            elif opt in ("-m", "--minlen"):
                minimumReadLength = int(arg)
            elif opt in ("-M", "--maxlen"):
                maximumReadLength = int(arg)
            
            elif opt in ("-q", "--minqual"):
                minimumQuality = int(arg)    
            
            elif opt in ("-r", "--rundate"):
                rundate = arg
            elif opt in ("-b", "--barcodes"):
                barcodeSampleMapFilename = arg
                barcodeSampleMap = readBarcodeSampleNames(barcodeSampleMapFilename)


    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    print('Input Directory:' + inputRootDirectory)
    #print('Output Directory:' + str(outputResultDirectory))
    print('Minimum Read Length:' + str(minimumReadLength))
    print('Maximum Read Length:' + str(maximumReadLength))

    # Quick sanity check.
    if(len(inputRootDirectory) < 4):
        print('Input directory is too short:' + str(inputRootDirectory))
        return False
    if(len(inputRootDirectory) < 4):
        print('Output directory is too short:' + str(outputResultDirectory))
        return False
    if(minimumReadLength < 0):
        print('Minimum Read Length should be >= 0 :' + str(minimumReadLength))
        return False
    
    if(minimumQuality < 0):
        print('Minimum Read Quality >= 0 :' + str(minimumQuality))
        return False
    # Actually maximumReadLength can be 0.  That's fine.  No max Length in that case  
    if(maximumReadLength < 0):
        print('Maximum Read Length should be >= 0 :' + str(maximumReadLength))
        return False
    if(rundate is not None and ('\\' in rundate or '/' in rundate or ' ' in rundate)):
        print('Please use normal characters in the rundate, no spaces or slashes. Underscore(_) and period(.) are okay.')
        return False
    
    #TODO: Check if there is a barcode file.  
    # If not, check for one in the root directory of where this is executing.
    # If not, make one in root dir. 

    return True

def usage():
    print("usage:\n" + 
    "\tThis script is written for python 2.7.11\n" + 
    "\tYou must specify both an input directory, and output directory to store the extracted reads.\n" + 
    "\t$python MinION_Read_Extractor_CL.py -i inputDirectory -o outputDirectory\n\n" + 
    "\tThe input directory may contain either (or both) of:\n" + 
    "\t\tAny number of .fast5 MinION reads (ex. a \"/pass\" folder)\n" + 
    "\t\tOR\n" + 
    "\t\tSubdirectories containing .fast5 reads (Barcoded Reads, ex. /BCO1)\n\n" + 
    "\tThe output directory will be filled with .fasta and .fastq reads, sorted by subfolders.\n\n" + 

    "\tOptions:\n" +  
    "\t-i\t--idir   \tInput Directory (required)\n" +  
    "\t-o\t--odir   \tOutput Directory (required)\n" +  
    "\t-m\t--minlen \tMinimum Read Length Filter\n" +  
    "\t-M\t--maxlen \tMaximum Read Length Filter\n" +  
    "\t-h\t--help   \tPrint this message\n" +   
    "\t-r\t--rundate\tSequencing Date or other information which will be included in the extract filename.\n" +   
    
    "\n\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
    )    

if __name__=='__main__':
    global inputRootDirectory
    global outputResultDirectory
    global minimumReadLength
    global maximumReadLength
    global minimumQuality
    global barcodeSampleMapFilename
    global barcodeSampleMap
    global rundate

    try:
    
        if(readArgs()):
            print('Commandline arguments look fine.  Starting to extract reads from:' + inputRootDirectory)
            
            myReadExtractor = MinION_Read_Extractor(inputRootDirectory, outputResultDirectory)
            
            myReadExtractor.minimumReadLength = minimumReadLength
            myReadExtractor.maximumReadLength = maximumReadLength        
            myReadExtractor.minimumQuality    = minimumQuality
            myReadExtractor.barcodeSampleMap  = barcodeSampleMap
            myReadExtractor.rundate           = rundate
            
            # Start extracting, and begin in the root directory.        
            myReadExtractor.extractAllFast5(myReadExtractor.inputRootDirectory, True)
            
            print ('Done extracting reads, have a nice day.')
    
        else:
            print('\nI give up, I\'m exiting without extracting your reads.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Fatal problem during read extraction:'
        print sys.exc_info()
        raise

