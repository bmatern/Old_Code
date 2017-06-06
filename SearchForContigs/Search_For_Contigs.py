import sys
import os
import subprocess
#from subprocess import check_output
from os import listdir, remove
from os.path import split, join, isfile, getsize
from shutil import copyfile, copy, copy2
#import random
#import math
from Bio import SeqIO
#from Bio import SeqRecord
#from Bio.Blast.Applications import NcbiblastnCommandline
#from Blast_Result import Blast_Result

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        os.mkdir(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

class Search_For_Contigs:
    inputDirectory = ''
    outputDirectory = ''
    contigResultsOutput = ''

    def WriteContigInfo(self, fileName, contigCount, contigLengths, contigReadCounts, contigCovStats):
        self.contigResultsOutput.write(
            str(fileName) + '\t' + 
            str(contigCount) + '\t' + 
            str(contigLengths) + '\t' + 
            str(contigReadCounts) + '\t' + 
            str(contigCovStats) + '\n')

    def SearchForContigs(self):
        print ('Looking for Contigs.')

        self.contigResultsOutput = createOutputFile(join(self.outputDirectory , 'Contigs_Found.txt')) 

        self.WriteContigInfo("File_Name","Contig_Count", "Contig_Lengths", "Contig_Read_Counts", "Contig_Cov_Stats")
        #self.contigResultsOutput.write('Contig files found:\n')

        # Start at the input directory, searching for contig files.
        self.SearchSubdirectory(self.inputDirectory)
           
        # Each record represents an HLA element in the input fasta file.
        #for index, record in enumerate(contigRecords):
        #    currentReadID = str(record.id)
            #print ('Read ID:' + currentReadID)
        #    currentSequence = str(record.seq)
            #print ('Read Sequence:' + currentSequence)

            #print ('Sorting Read (' + str(index) + '/' + str(readCount) + ') : ' + currentReadID)

        self.contigResultsOutput.close()

    def SearchSubdirectory(self, subDirectory):
        #self.contigResultsOutput.write('Searching Subdirectory:' + subDirectory + '\n')

        fileNames = [f for f in listdir(subDirectory) if isfile(join(subDirectory, f))]
        fileNames = sorted(fileNames)
        #self.contigResultsOutput.write('FileCount:' + str(len(fileNames)) + '\n')
        for fileName in fileNames:
            if('contigs.fasta' in fileName):
                #self.contigResultsOutput.write('Contig File Found: ' + fileName + '\n')
                self.FoundContig(fileName, subDirectory)             
            elif('bubbles.fasta' in fileName): 
                #self.contigResultsOutput.write('Bubble File Found: ' + fileName + '\n')
                self.FoundContig(fileName, subDirectory)
            else:
                pass

        folderNames = [f for f in listdir(subDirectory) if not isfile(join(subDirectory, f))]
        folderNames = sorted(folderNames)
        #self.contigResultsOutput.write('FolderCount:' + str(len(folderNames)) + '\n')
        for folderName in folderNames:
            #self.contigResultsOutput.write('Recursively calling for subdir:' + join(subDirectory,folderName) + '\n')
            self.SearchSubdirectory(join(subDirectory,folderName))

        pass

    def FoundContig(self, fileName, subDirectory):
        fileSize = getsize(join(subDirectory,fileName))
        #self.contigResultsOutput.write('FileSize:' + str(fileSize) + '\n')
        
        if(fileSize > 0):

            #self.contigResultsOutput.write('Copying this file to output directory\n')
            if isfile(join(self.outputDirectory,fileName)):
                #File already exists.
                pass
            else:
                copy2( join(subDirectory,fileName) , self.outputDirectory)

            #Load Fasta file.  Count contigs.
            # Each record represents an HLA element in the input fasta file.

            contigLengths = '{'
            contigReadCounts = '{'
            contigCovStats = '{'

            contigRecords = SeqIO.parse(join(subDirectory,fileName), "fasta")
            contigCount = len(list(SeqIO.parse(join(subDirectory,fileName), "fasta")))
            for index, record in enumerate(contigRecords):
                #self.contigResultsOutput.write('Record:' + str(record) + '\n')

                currentDescription = str(record.description)
                # Description looks like this:
                # tig00000000 len=3161 reads=59 covStat=312.07 gappedBases=no class=contig suggestRepeat=no suggestCircular=no
                headerTokens = currentDescription.split()

                #lengthToken = headerTokens[1]
                currentContigLength = headerTokens[1][headerTokens[1].find('=') + 1 : len(headerTokens[1])]
                contigLengths += currentContigLength + ';'
                #self.contigResultsOutput.write('Contig Length:' + currentContigLength + '\n')

                #readCountToken = headerTokens[2]
                currentContigReadCount = headerTokens[2][headerTokens[2].find('=') + 1 : len(headerTokens[2])]
                contigReadCounts += currentContigReadCount + ';'
                #self.contigResultsOutput.write('Contig Read Count:' + currentContigReadCount + '\n')

                currentContigCovStat = headerTokens[3][headerTokens[3].find('=') + 1 : len(headerTokens[3])]
                contigCovStats += currentContigCovStat + ';'
                #self.contigResultsOutput.write('Contig Cov Stat:' + currentContigCovStat + '\n')

            contigLengths = contigLengths[0:len(contigLengths)-1] + '}'
            contigReadCounts = contigReadCounts[0:len(contigReadCounts)-1] + '}'
            contigCovStats = contigCovStats[0:len(contigCovStats)-1] + '}'

            self.WriteContigInfo(fileName,contigCount,contigLengths,contigReadCounts,contigCovStats)
             

        else:
            #self.contigResultsOutput.write('It\'s empty, skipping this file.\n')
            self.WriteContigInfo(fileName,"0", "-", "-", "-")
            pass


if __name__=='__main__':
    try:
        ContigSearcher = Search_For_Contigs()

        ContigSearcher.inputDirectory = sys.argv[1]
        ContigSearcher.outputDirectory = sys.argv[2]
        ContigSearcher.SearchForContigs()
        print('Done.  Yay.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise







