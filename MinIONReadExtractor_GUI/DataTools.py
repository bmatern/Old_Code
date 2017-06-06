# This file is part of MinION-extractor-GUI.
#
# MinION-extractor-GUI is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MinION-extractor-GUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MinION-extractor-GUI. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0 

SoftwareVersion = "MinION-extractor-GUI_Version_1.0"

import h5py
import numpy
from Bio import SeqIO
from StringIO import StringIO
import os
from os import listdir
from os.path import isfile, join, split
import tkMessageBox
from MinIONRead import MinIONRead

# This is the location of the two-directional data within a fast5 file.  
# This might change, depending on Oxford Nanopore file format changes.
# Perhaps it's possible to extract additional information by messing with this/these keys.  
keys = {'twodirections' : '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'}


# The MinIONReadWriter class contains logic for extracting sequences from fast5 files
# and writing the sequences to fasta or fastq output.
class MinIONReadWriter:

    # This method is a directory-safe way to open up a write file.
    @staticmethod
    def createOutputFile(outputfileName):
        tempDir, tempFilename = os.path.split(outputfileName)
        if not os.path.isdir(tempDir):
            os.mkdir(tempDir)
        resultsOutput = open(outputfileName, 'w')
        return resultsOutput

    # writeResultsOutput method calculates simple stats on the reads, and writes them to a file.  
    @staticmethod
    def writeResultsOutput(reads, readResultsFileName):

        #Lets skip this step if there are no reads.
        if(len(reads) > 0):

            resultsOutput = MinIONReadWriter.createOutputFile(readResultsFileName)
            reads.sort(key=lambda x: x.sequenceLength, reverse=False)
        
            resultsOutput.write('ReadCount:' + str(len(reads)) 
            + ' Min:' + str(reads[0].sequenceLength)        
            + ' Max:' + str(reads[len(reads) - 1].sequenceLength)
            + ' Mean Length:' + str(sum(read.sequenceLength for read in reads) /len(reads))
            + ' Mean Quality:' + str(sum(read.averageQuality for read in reads) /len(reads))
            + '\n\n')
        
            resultsOutput.write('Length\tQuality\n')

            for i in range(0,len(reads)):
                resultsOutput.write(str(reads[i].sequenceLength) 
                + '\t' + str(reads[i].averageQuality)
                + '\n')

            resultsOutput.close()


    # writeOutput contains logic to read fast5 files, and write the reads to output fasta and fastq.
    # This method returns a boolean
    # Which represents the answer to the question "Should we destroy this window and continue?"
    # Should almost always be true.
    @staticmethod
    def writeOutput(inputDirectory, writeFasta, fastaFileName, writeFastq, fastqFileName, writeReadResults, readResultsFileName):
        
        if (len(inputDirectory) < 1):
            tkMessageBox.showinfo(title='No Input Directory', message='Select an input directory containing MinION .fast5 files.')

        elif (not writeFastq and not writeFasta):
            tkMessageBox.showinfo(title='No Output File', message='You must output to a .fastq or .fasta file.')

        else:

            # Ask confirmation before starting
            confirmationMessage = ('You will read .fast5 files from:\n'
            + inputDirectory + '\n\n')
            if  (writeFastq):
                confirmationMessage = (confirmationMessage + 'Write to this .fastq file:\n'
                + fastqFileName + '\n\n')
            if  (writeFasta):
                confirmationMessage = (confirmationMessage + 'Write to this .fasta file:\n'
                + fastaFileName + '\n\n')   
            if  (writeReadResults):
                confirmationMessage = (confirmationMessage + 'Write read stats to:\n'
                + readResultsFileName + '\n')       
     
            if (tkMessageBox.askokcancel(title='Are you sure?', message=confirmationMessage)):

                fileNames = [f for f in listdir(inputDirectory) if isfile(join(inputDirectory, f))]
                if(len(fileNames) > 0):
                    if  (writeFastq):
                        outputFileFastq = MinIONReadWriter.createOutputFile(fastqFileName)

                    if  (writeFasta):
                        outputFileFasta = MinIONReadWriter.createOutputFile(fastaFileName)
            
                    currentSequence = 0
                    totalCount = len(fileNames)
                    
                    msg = ('Extracting sequences from ' + str(totalCount) + ' .fast5 files in directory:\n' 
                    +  inputDirectory + '\nPlease wait...')
                    print(msg)

                    # List of simple read stats
                    minIONReadData = []
                    
                    # Loop through each file in this directory.
                    for fileName in fileNames:
                        try:
                            currentSequence+=1
                            
                            combinedFileName = join(inputDirectory, fileName)

                            #Progress text.  This is printed to console.
                            if( (currentSequence == 1) or (currentSequence == totalCount) or (currentSequence % 500)==0 ):
                                percentage = str( (100.00 * currentSequence) / totalCount )[:6] + '%'
                                msg = ('Reading file: '
                                + '(' + str(currentSequence) + '/' + str(totalCount) + ') ... ' + percentage)
                                print(msg)	
                         
                            # Only want fast5 files.  
                            if('fast5' in fileName):		
                                hdf = h5py.File(combinedFileName, 'r')
                                for id, key in keys.iteritems():    

                                    fq = hdf[key][()]

                                    rec = SeqIO.read(StringIO(fq), "fastq")
                                    rec.id += "_" + id + " " + SoftwareVersion

                                    # "qualities" will contain letter-by-letter quality scores 
                                    # for each sequenced nucleotides. 
                                    # Mankind has not yet discovered how ONT calculates quality scores.
                                    # An algebraic mean doesn't makes much sense for logarithmic phred scores.
                                    # It's something.  
                                    qualities = rec.letter_annotations["phred_quality"]
                                    currentAvgQuality = numpy.mean(qualities)

                                    currentMinIONRead = MinIONRead()
                                                                    
                                    currentMinIONRead.id = rec.id
                                    currentMinIONRead.sequenceLength = len(rec)
                                    currentMinIONRead.averageQuality = currentAvgQuality
                                    
                                    minIONReadData.append(currentMinIONRead)

                                    if (writeFastq):
                                        SeqIO.write([rec], outputFileFastq, 'fastq')
                                    if (writeFasta):
                                        SeqIO.write([rec], outputFileFasta, 'fasta')
                                hdf.close()
                        except Exception, e:
                            print('Exception occured when reading a .fast5 file: ' + str(e))
                            print('Skipping this file: ' + fileName)
                    
                    if (writeFastq):    
                        outputFileFastq.close()
                    if (writeFasta):
                        outputFileFasta.close()

                    if  (writeReadResults):
                        MinIONReadWriter.writeResultsOutput(minIONReadData, readResultsFileName)
                    
                # There are no files found in this directory.
                else:
                    tkMessageBox.showinfo(title='No files in this directory', 
                        message='There are no files in the folder\n'
                        + inputDirectory                        
                        + '\nSkipping this folder.')


                return tkMessageBox.askyesno(title='Done.  Close and continue now?', message=(
                    'Done extracting this directory.  ' +
                    'Do you want to close this window and continue now?'))

      
