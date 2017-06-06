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

import h5py
import numpy
from numpy import asarray

from Bio import SeqIO
from StringIO import StringIO
import os
from os import listdir
from os.path import isfile, join, split
import sys
import traceback
#import subprocess
#import MinIONRead
import MinIONReadOutputFile
import pylab
from matplotlib import collections

# This script can be used to extract the fastq and fasta sequences 
# from a directory of fast5 files, from MinION Output
# It expects a single parameter, to specify a search directory.  
# Output files are generated in a read_extracts folder under the input directory.
class MinION_Read_Extractor:
    
    def __init__(self, inputRootDirectory=None, outputResultDirectory=None):
        # These keys represent the locations of fastq data within a fast5 file.  
        # This might change, depending on Oxford Nanopore file format changes.
        self.fastqKeys = {
             '2D'           : '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'
            ,'1D_Template'   : '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
            ,'1D_Complement' : '/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq'
            ,'Barcoding'        : '/Analyses/Barcoding_000/Barcoding/Fastq'
        }
        
        # These are the locations of the metadata attributes on the Fast5 files.
        # Keep this in sync with metadataAttributeNames
        self.metadataLocations={
             'BarcodingChimaeraVersion'     : '/Analyses/Barcoding_000'
            ,'BarcodingDragonetVersion'     : '/Analyses/Barcoding_000'
            ,'Basecall1DChimaeraVersion'    : '/Analyses/Basecall_1D_000'
            ,'Basecall1DDragonetVersion'    : '/Analyses/Basecall_1D_000'
            ,'Basecall2DChimaeraVersion'    : '/Analyses/Basecall_2D_000'
            ,'Basecall2DDragonetVersion'    : '/Analyses/Basecall_2D_000'
            ,'CalibrationChimaeraVersion'   : '/Analyses/Calibration_Strand_000'
            ,'CalibrationDragonetVersion'   : '/Analyses/Calibration_Strand_000'
            ,'EventDetectionMinKNOWVersion' : '/Analyses/Hairpin_Split_000'
            ,'HairpinSplitChimaeraVersion'  : '/Analyses/Hairpin_Split_000'
            ,'HairpinSplitDragonetVersion'  : '/Analyses/EventDetection_000'
            ,'DeviceID'                     : '/UniqueGlobalKey/tracking_id'
            ,'ProtocolsVersionName'         : '/UniqueGlobalKey/tracking_id'
            ,'TrackingVersion'              : '/UniqueGlobalKey/tracking_id'
        }
    
        # These are the names of the attributes of the metadata\
        # keep this in sync with metadataLocations
        self.metadataAttributeNames={
            'BarcodingChimaeraVersion'      : 'chimaera version'
            ,'BarcodingDragonetVersion'     : 'dragonet version'
            ,'Basecall1DChimaeraVersion'    : 'chimaera version'
            ,'Basecall1DDragonetVersion'    : 'dragonet version'
            ,'Basecall2DChimaeraVersion'    : 'chimaera version'
            ,'Basecall2DDragonetVersion'    : 'dragonet version'
            ,'CalibrationChimaeraVersion'   : 'chimaera version'
            ,'CalibrationDragonetVersion'   : 'dragonet version'
            ,'HairpinSplitChimaeraVersion'  : 'chimaera version'
            ,'HairpinSplitDragonetVersion'  : 'dragonet version'
            ,'EventDetectionMinKNOWVersion' : 'version'
            ,'DeviceID'                     : 'device_id'
            ,'ProtocolsVersionName'         : 'protocols_version_name'
            ,'TrackingVersion'              : 'version'
        }
        
        # Display progress text every x reads.  250 is a good choice?
        self.progressBarBlockSize=250

        self.inputRootDirectory    = inputRootDirectory
        self.outputResultDirectory = outputResultDirectory
        self.minimumReadLength     = 0
        self.maximumReadLength     = 0
        self.minimumQuality        = 0
        self.barcodeSampleMap      = None
        self.rundate               = None
    
        self.outputFastaFiles      = []
        self.outputFastqFiles      = []
        
        self.rejectLengthOutput    = None
        self.rejectQualityOutput   = None
        self.readListOutput        = None
        
        self.readListOutputFileName= None
        
        # A small place to keep track of metadata.
        # A dictionary so it's easy to look things up.
        self.metadataAttributes    = {}

        
    def createScatterPlot(self, graphTitleText, xValues, yValues, xAxisName, yAxisName, outputFileName):
        print('Creating a Scatter Plot: ' + outputFileName)
        
        #Clear the figure and start anew.
        pylab.clf()
        
        # K is black, we need to repeat N times.
        colors = ['K'] * len(xValues)
    
        pylab.scatter(xValues, yValues, s=1, c=colors, marker='.')
        
        pylab.xlabel(xAxisName)
        pylab.ylabel(yAxisName)
        pylab.title(graphTitleText)
        
        pylab.savefig(join(self.outputResultDirectory, outputFileName))
           

    # Short method to print out summary stats for a group of reads
    def writeReadStatsSub(self, readType, barcode, readStats):
        
        statsSubText=''        
        
        if(len(readStats) > 0):
            #readLengths = readStats[:,0]
            #readQualities = readStats[:,1]
            
            statsSubText += ('\n' + barcode + ' ' + readType + ' Reads Summary:\n' 
                + 'Sequences Extracted,' + str(len(readStats)) + '\n'
                + 'Minimum Length,' + str(int(numpy.amin(readStats[:,0]))) + '\n'
                + 'Maximum Length,' + str(int(numpy.amax(readStats[:,0]))) + '\n'
                + 'Mean Length,' + str(numpy.mean(readStats[:,0])) + '\n'
                + 'Mean Quality,' + str(numpy.mean(readStats[:,1])) + '\n'
            )
            
        else:
            statsSubText += ('\n' + barcode + ' ' + readType + ' Reads Summary:\n' 
                + 'Sequences Extracted,0\n'

            )

        return statsSubText
    
    # Recursive method to get a list of all fast5 files in a directory
    def findAvailableReads(self, currentDirectory, recurseSubdirectories):
        print 'Looking for fast5 reads ... ' + currentDirectory
        
        # Get children files and folders from the input currentDirectory
        # This will only look at .fast5 files.
        readFileNames = [f for f in listdir(currentDirectory) if (isfile(join(currentDirectory, f)) and f.endswith('.fast5'))]
        readFileNames.sort()
        
        subdirectoryNames = [f for f in listdir(currentDirectory) if not isfile(join(currentDirectory, f))]
        subdirectoryNames.sort()
        
        # Use a set, so it will automatically get rid of duplicates.
        readsToExtract = set()
        
        # Add all the files in *this* directory.
        for readFileName in readFileNames:
            readFullPath = os.path.join(currentDirectory, readFileName)
            
            #print('Looking for reads, foundthisone:\n' + readFullPath + '\n*')
            readsToExtract.add(readFullPath)
                 
        #Recurse into each subdirectory.
        # Add all the files from each subdirectory.
        if(recurseSubdirectories):
            for subdirectoryName in subdirectoryNames:
                    
                subdirFullPath = os.path.join(currentDirectory, subdirectoryName)
                subdirReads = self.findAvailableReads(subdirFullPath, recurseSubdirectories)
                
                readsToExtract = readsToExtract.union(subdirReads)
  
        return readsToExtract
       
    # This constructs a list of the .fast5 reads that we have ALREADY extracted.             
    
    
    def getReadsAlreadyExtracted(self):        
        readList=[]
        try:
            print 'Checking the output folder extracted reads : ' + self.outputResultDirectory
            
            if self.readListOutputFileName is not None:
            
                print('Try to close the file list')
                try:
                    self.readListOutput.close()
                except Exception, e:
                    print('Cannot close it, maybe already closed, maybe doesn\'t exist. Either way no problem.')
              
                if  os.path.exists(self.readListOutputFileName):            
                    print 'The Read Output File is here: ' + self.readListOutputFileName
                    
                    with open(self.readListOutputFileName, 'r') as readListReader:
                        for line in readListReader.readlines():
                            readList.append(line.strip())
                        #readList = readListReader.readlines()
                        
                    #here im reading the output file.  Should check if i tis open first.  
    
                    return readList  
                else:
                    print('Extracted Read File doesn\'t exist, i\'ll have to return an empty set.')
                    return readList
            else:
                print ('Read Output File does not exist yet, so I won\'t read it.')
                return readList
        
        except Exception, e:
            print('Exception when reading the list of already extracted reads:  ' + str(e))
            # TODO: I really should check if the file is open already, 
            # I think that could throw an exception here.
            raise

                    
    def extractAllFast5(self, rootDirectory, recurseSubdirectories):
        print('Step 0: Validate all my input parameters.  are all the directoreies set? I still need to do this.')
        # TODO : Validate all my input parameters.  Do i have a barcode map? 
        # Do I have an input directory?  
        # Do I have an ouptput directory? Does it exist?  Create it.  This is necessary.
        
        if not os.path.exists(self.outputResultDirectory):
            print('Creating the Output Directory: ' + self.outputResultDirectory)
            os.makedirs(self.outputResultDirectory)
        else:
            print('Not Creating the Output Directory because it exists already: ' + self.outputResultDirectory)
      
        print('Step 0.5: Create Output Files')
        self.createSequenceFileOutputs()
            
        
        print('Step 1: Get a list of all the reads we can extract.')  
        availableReads = self.findAvailableReads(rootDirectory, recurseSubdirectories)
       
        print('Step 2: What reads are already extracted?')  
        readsAlreadyExtracted = self.getReadsAlreadyExtracted()
      
        print('Step 3: Skip the reads we have already extracted.')
        print(str(len(availableReads)) + ' reads available')
        #print('for example this read:' + str(availableReads.pop()))
        print(str(len(readsAlreadyExtracted)) + ' reads already completed')
        #print('for example this read:' + str(readsAlreadyExtracted.pop()))
        
        newReads = [read for read in availableReads if read not in readsAlreadyExtracted]
        print(str(len(newReads)) + ' reads ready for extraction now')
        
        print('Step 4: Process new reads, OR, finish up and draw graphs')
        # If we have no new reads to process.
        if(len(newReads) == 0):
            print('No new reads for extracting were found.')
            self.calculateSummariesAndFinish()
            # Make Graphs
            # Make Summary Files, by opening the output files and summarizing.
            pass
        else:
            # Extract the new reads
            
            # Then, recurse because there might be NEW reads to process(If MinION is still running, for example)
            if(self.extractNewReads(newReads)):
                print('Done extracting new reads.  I will now check for new reads to make sure I got them all.')
                # TODO: Sleep for 5 seconds or something before trying again.
                # TODO: Maybe I should generate graphs quick too, to make it realtime.  Probably not though.
                # self.calculateSummaries()  This method might not exist yet.
                # Depends on how quick it is.
                
                # TODO: Maybe there should be a maximum duration parameter.  Extract for maximum 1 hour?
                self.extractAllFast5(rootDirectory, recurseSubdirectories)
              
            else:
                print('There was a problem when extracting the new reads.  I will not try again, please debug this.')
            pass
        
        print('Done extracting fast5 from root directory : ' + rootDirectory)
        
        
    # This method will calculate an array of read Lengths and Qualities
    # From a fastq filename
    def getReadStats(self, inputFastqFileName): 
        try:
            print('Calculating read stats for this file:' + str(inputFastqFileName))
            parsedExtractedReads = SeqIO.parse(inputFastqFileName, 'fastq')
            minionReadRecords = enumerate(parsedExtractedReads)
            
            sequenceStats = numpy.empty((0,2), float)
            
            for index, record in minionReadRecords: 
    
                qualities = record.letter_annotations["phred_quality"]                    
                currentSeqLength = len(record)
                currentAvgQuality = numpy.mean(qualities)
                
                sequenceStats = numpy.append(sequenceStats, numpy.array([[currentSeqLength, currentAvgQuality]]), axis=0 )

            return sequenceStats
        
        except Exception, e:
            print 'Exception trying to calculate read stats for this file:' + str(inputFastqFileName)
            print sys.exc_info()
            raise
        
    def calculateSummariesAndFinish(self):
        print('Summarizing Extracted Reads and Generating Read Scatterplots.')

        # TODO Open the output result file here.  Don't give it a barcode name.
        # Open for writing, so it overwrites whatever is in there before.  
        summaryFileName = os.path.join(self.outputResultDirectory,  'ReadSummary.csv') 
        summaryOutputFile = open(summaryFileName,'w')
        
        summaryOutputFile.write('Summary of .fast5 Read Extraction\n\n')
        
        summaryOutputFile.write('Extraction Parameters:\n')
        summaryOutputFile.write('Input Directory,' + self.inputRootDirectory + '\n')
        summaryOutputFile.write('Result Directory,' + self.outputResultDirectory + '\n')
        summaryOutputFile.write('Run Date,' + str(self.rundate) + '\n')
        summaryOutputFile.write('Minimum Length,' + str(self.minimumReadLength) + '\n')
        summaryOutputFile.write('Maximum Length,' + str(self.maximumReadLength) + '\n')
        summaryOutputFile.write('Minimum Quality,' + str(self.minimumQuality) + '\n\n')

        print('This many barcodes were found to summarize, including unbarcoded:' + str(len(self.outputFastqFiles)))
        
        # This is a text buffer so we can print things in the correct order. 
        barcodeSummaryTextBuffer=''
        
        allExtractedReadStats = None
        
        # for each barcode type in my output file 2d array
        for barcodeIndex in range(0,len(self.outputFastqFiles)):
            
            barcodeText = self.outputFastqFiles[barcodeIndex][0].barcode
            
            # foreach readtype
            for readTypeIndex in range(0,len(self.outputFastqFiles[barcodeIndex])):
                fastqFileName = self.outputFastqFiles[barcodeIndex][readTypeIndex].fileName 
                scatterplotFileName = fastqFileName.replace('.fastq','.png')
                
                readTypes = self.fastqKeys.keys()
                readTypes.sort()
                readType = str(readTypes[readTypeIndex])
                scatterplotTitle = 'Q vs. Length: ' + barcodeText + ' ' + readType + ' reads' 
                

                print('Summarizing, Barcode=' + barcodeText + ' , ReadType=' + str(readTypeIndex))
                print('This data is contained in the fastq file:' + fastqFileName)
                print('I\'ll make a scatterplot:' + scatterplotFileName)
                
                fastqReadStats = self.getReadStats(fastqFileName)
          
                if(allExtractedReadStats is None or len(allExtractedReadStats) == 0):
                    allExtractedReadStats = fastqReadStats
                else:
                    allExtractedReadStats = numpy.vstack((allExtractedReadStats,fastqReadStats))

                #self.createScatterPlot(scatterplotTitle, sequenceLengths, sequenceQualities, 'Sequence Length (bp)', 'Sequence Quality (phred)', scatterplotFileName)
                barcodeSummaryTextBuffer += self.writeReadStatsSub(readType, barcodeText, fastqReadStats)
                
                if(len(fastqReadStats)>0):
                    self.createScatterPlot(scatterplotTitle, fastqReadStats[:,0], fastqReadStats[:,1], 'Sequence Length (bp)', 'Sequence Quality (phred)', scatterplotFileName)
             
             
        # Header
        summaryOutputFile.write('Overall Summary:\n')
        summaryOutputFile.write('Type,Count,AverageLength,AverageQuality\n')
             
        # Calculate Length of the read list.  This should be the number of .fast5 files analyzed
        lineCount = 0
        try:
            with open(self.readListOutputFileName, 'r') as readListReader:
                lineCount = len(readListReader.readlines())
            readListReader.close()
            summaryOutputFile.write('.fast5 Reads Analyzed,' + str(lineCount) + '\n')            
        except Exception, e:
            print 'Exception trying to count the reads we extracted:'
            print sys.exc_info()
            raise
            
        lengthRejectedReadStats = self.getReadStats(self.rejectLengthOutput.fileName)
        qualityRejectedReadStats = self.getReadStats(self.rejectQualityOutput.fileName)
        
        rejectedReadStats = numpy.vstack((qualityRejectedReadStats,lengthRejectedReadStats))        
        allCombinedSequenceStats = numpy.vstack((rejectedReadStats,allExtractedReadStats))
        
        if(len(allCombinedSequenceStats) > 0):
            summaryOutputFile.write('Total Sequences Found,' + str(len(allCombinedSequenceStats)) + ',' + str(numpy.average(allCombinedSequenceStats[:,0])) + ',' + str(numpy.average(allCombinedSequenceStats[:,1])) + '\n')
        else:
            summaryOutputFile.write('Total Sequences Found,0\n')
        
        if(len(allExtractedReadStats) > 0):
            summaryOutputFile.write('Sequences Extracted,' + str(len(allExtractedReadStats)) + ',' + str(numpy.average(allExtractedReadStats[:,0])) + ',' + str(numpy.average(allExtractedReadStats[:,1])) + '\n')
        else:
            summaryOutputFile.write('Sequences Extracted,0\n')
        
        if(len(lengthRejectedReadStats) > 0):
            summaryOutputFile.write('Sequences Rejected for Length,' + str(len(lengthRejectedReadStats)) + ',' + str(numpy.average(lengthRejectedReadStats[:,0])) + ',' + str(numpy.average(lengthRejectedReadStats[:,1])) + '\n')
        else:
            summaryOutputFile.write('Sequences Rejected for Length,0\n')
            
        if(len(qualityRejectedReadStats) > 0):
            summaryOutputFile.write('Sequences Rejected for Quality,' + str(len(qualityRejectedReadStats)) + ',' + str(numpy.average(qualityRejectedReadStats[:,0])) + ',' + str(numpy.average(qualityRejectedReadStats[:,1])) + '\n')
        else:
            summaryOutputFile.write('Sequences Rejected for Quality,0\n')
        
         
        # barcodeSummaryTextBuffer contains information about the extracted read types.
        summaryOutputFile.write(barcodeSummaryTextBuffer)
        
        summaryOutputFile.write(self.writeReadStatsSub('Length Rejected', 'All Barcodes', lengthRejectedReadStats))
        summaryOutputFile.write(self.writeReadStatsSub('Quality Rejected', 'All Barcodes', qualityRejectedReadStats))

        
        if(len(lengthRejectedReadStats)>0):
            self.createScatterPlot('Q vs. Length, Length-Rejected reads' , 
                lengthRejectedReadStats[:,0], 
                lengthRejectedReadStats[:,1], 
                'Sequence Length (bp)', 'Sequence Quality (phred)'
                , self.rejectLengthOutput.fileName.replace('.fastq','.png'))
            
        if(len(qualityRejectedReadStats)>0):
            self.createScatterPlot('Q vs. Length, Quality-Rejected reads' , 
                qualityRejectedReadStats[:,0], 
                qualityRejectedReadStats[:,1], 
                'Sequence Length (bp)', 'Sequence Quality (phred)'
                , self.rejectQualityOutput.fileName.replace('.fastq','.png'))

             
        # Write metadata          
        for metaDataKeyId, metaDataKey in self.metadataLocations.iteritems():        
            try:        
                summaryOutputFile.write('\nRead Metadata:\n')

            except Exception:
                # Top Level exception handling like a pro.
                # This is not really doing anything.
                print 'Exception trying to write metadata to results file:'
                print sys.exc_info()
                traceback.print_exc()
                
        # Close output file.   

        summaryOutputFile.close()

    # This method will create output files for a specific Barcode found within a MinION read.
    # The return value is the index of the barcode within the outputFile2dList
    def addOutputBarcode(self, outputFile2dList, barcodeText, fileFormat):
        # If barcode exists
        for barcodeIndex in range(0,len(outputFile2dList)):
            if outputFile2dList[barcodeIndex][0].barcode == barcodeText:
                # We found the barcode already.  No need to continue.
                return barcodeIndex
            
        print('Adding Barcode ' + str(barcodeText) + ' to the ' + str(fileFormat) + ' output files.')
        
        # If we escaped that loop, then this barcode's output files
        # don't exist yet.  Create, Add, and Open the output files.
        # Open them for appending, not writing.
        readTypes = self.fastqKeys.keys()
        readTypes.sort()
        newBarcodeRow = [None]*len(readTypes)

        for readTypeIndex, readType in enumerate(readTypes):
            newBarcodeRow[readTypeIndex] = MinIONReadOutputFile.MinIONReadOutputFile(self.rundate, barcodeText, readType, self.outputResultDirectory, fileFormat)
        
        outputFile2dList.append(newBarcodeRow) 
        
        return len(outputFile2dList) - 1

        
    

    def extractSingleMinIONRead(self, readFilename):      
        #print 'Extracting Read : ' + str(readFilename)
        # Open .fast5 file.
        try:
            hdf = h5py.File(readFilename, 'r')      
            
            # foreach read type  
            # Sort the read types to stay consistent   
            readTypes = self.fastqKeys.keys()
            readTypes.sort() 
            
            #print('entering loop')
            for readTypeIndex, readType in enumerate(readTypes):
            
                #print('Fast5 Open, key id:' + id)
            
                #sequenceFound = False
                rec = None
                # Attempt to read the fastq file at this location
                try:
                    fq = hdf[self.fastqKeys[readType]][()]
                    
                    #print('fast5 is opened.')
            
                    rec = SeqIO.read(StringIO(fq), "fastq")                    
                    rec.id += "_" + readType + " " + SoftwareVersion
                    # Add the location of the read file in the header.  
                    # nanopolish expects to find the file location in the last token of the fasta header.
                    #https://github.com/jts/nanopolish/blob/45335166365ccc46a3dd434f5b937f1594e103f8/src/common/nanopolish_fast5_map.cpp
                    rec.description += " " + readFilename
                    
                    qualities = rec.letter_annotations["phred_quality"]
                    currentAvgQuality = numpy.mean(qualities)
                    
                    #print('checking lengths')
                    
                    # Fetch Metadata.
                    #self.fetchMetaData(hdf)
                    # TODO: Uncomment the call for fetchMetaData.  The method needs work.    
                    
                    if( (len(rec.seq) < self.minimumReadLength and self.minimumReadLength != 0)
                        or (len(rec.seq) > self.maximumReadLength and self.maximumReadLength != 0) 
                        ):  
                        #print('Read length (' + str(len(rec.seq)) + ') is outside specified parameters. I will reject this read for Length.')   
                        SeqIO.write([rec], self.rejectLengthOutput, 'fastq')
                        
                        
                    # If read is below quality
                    elif( currentAvgQuality < self.minimumQuality and self.minimumQuality != 0
                        ):          
                        #print('Read quality is too low.  I will reject this read for quality.')
                        SeqIO.write([rec], self.rejectQualityOutput, 'fastq')
                    
                    # Else we have a good quality read that we want to extract.
                    else:
                        # Choose output file from 2d array
                        
                        # Get Barcode.  From the .fast5 file, not from the directory name, foolish mortal.
                        # The barcode lives in the fast5 file here:
                        # '/Analyses/Barcoding_000/Barcoding/barcode_arrangement'
                        # TODO: Barcode location has moved.  I see it here now:
                        # /Analyses/Barcoding_000/Summary/barcoding
                        try:
                            barcodeKey = '/Analyses/Barcoding_000/Barcoding'
                            currentBarcode = hdf[barcodeKey].attrs['barcode_arrangement']
                
                            #currentBarcode=hdf[barcodeKey]#[()]
                        except KeyError:
                            #print('barcode key exception...')
                            currentBarcode='unbarcoded'
                            
                        # Get barcode Index
                        # This will check for duplicates, no worries.
                        barcodeOutputIndex = self.addOutputBarcode(self.outputFastaFiles, currentBarcode, 'fasta')
                        barcodeOutputIndex = self.addOutputBarcode(self.outputFastqFiles, currentBarcode, 'fastq')         
                        
                        # Append my read to the end of the file.
                        currentFastaOutputFile = self.outputFastaFiles[barcodeOutputIndex][readTypeIndex].outputFile
                        currentFastqOutputFile = self.outputFastqFiles[barcodeOutputIndex][readTypeIndex].outputFile
                        
                        SeqIO.write([rec], currentFastaOutputFile, 'fasta')
                        SeqIO.write([rec], currentFastqOutputFile, 'fastq')
                    
                    
            
                except KeyError:
                    # These exceptions are completely expected and okay.  
                    # There is no data in this key, and that's fine most of the time.
                    # I am masking potential problems by ignoring exceptions, which makes me nervous.
                    pass

            hdf.close()
        except Exception, e:
            print('Exception when opening .fast5 file: ' + str(readFilename) + '\nProbably ok if this doesn\'t happen often: ' + str(e))
            print('If you see this message alot, you should look into what is wrong with all these files.')
            
            # TODO: Sometimes fast5 files are corrupted.  I should handle this case a little better. 
            #raise()

    def createSequenceFileOutputs(self):
        # Create output files.
        # 2x 2D list, for fasta and fastq
        
        # Sanity check, because if they already exist, I don't want to do this.        
        if(self.outputFastqFiles is None or len(self.outputFastqFiles) < 1):
            print('Creating Sequence File Outputs.')
    
            self.outputFastaFiles = []
            self.outputFastqFiles = []
            
            # The entries at index zero is for unbarcoded reads
            self.addOutputBarcode(self.outputFastaFiles, 'unbarcoded', 'fasta')
            self.addOutputBarcode(self.outputFastqFiles, 'unbarcoded', 'fastq')

        
        
        if(self.rejectLengthOutput is None):
            print('Creating the Length Reject Output file.')
            self.rejectLengthOutput = MinIONReadOutputFile.MinIONReadOutputFile(self.rundate, 'length_rejected', 'allreadtypes', self.outputResultDirectory, 'fastq')
            #print 'This is the length reject output file:' + str(self.rejectLengthOutput)
        else:
            print('Not creating the Length Reject Output file because it already exists.')
            
        if(self.rejectQualityOutput is None):
            print('Creating the Quality Reject Output file.')
            self.rejectQualityOutput = MinIONReadOutputFile.MinIONReadOutputFile(self.rundate, 'quality_rejected', 'allreadtypes', self.outputResultDirectory, 'fastq')
        else:
            print('Not creating the Quality Reject Output file because it already exists.')
                
        #Create the read list output file.       
        # I need to check if this object is closed.  
        if(self.readListOutput is None or self.readListOutput.closed):
            print('Creating the Read List Output file.')

            shortFileName = str(self.rundate) + '_' + 'extracted_reads.txt' 
            self.readListOutputFileName = os.path.join(self.outputResultDirectory, shortFileName)
            self.readListOutput = open(self.readListOutputFileName ,'a') 
            
            #print 'This is the length reject output file:' + str(self.rejectLengthOutput)
        else:
            print('Not creating the Read List Output file because it already exists.')
      
    def closeOutputFiles(self):
        # Close the 2x 2d arrays of output files.
        # I hope that the fasta and fastq files are still in sync, that makes this easier.

        print 'Closing the output read files.'
        if (len(self.outputFastaFiles) > 0):
            for barcodeIndex in range(0, len(self.outputFastaFiles)):
                for readTypeIndex in range(0, len(self.outputFastaFiles[0])):
                    self.outputFastaFiles[barcodeIndex][readTypeIndex].close()
                    self.outputFastqFiles[barcodeIndex][readTypeIndex].close()

        
        print ('Closing the Rejected Read outputs.')
        self.rejectLengthOutput.close()
        self.rejectQualityOutput.close()
        
        self.readListOutput.close()


    def extractNewReads(self, readList):
        try:
            readList.sort()
            print ('Extracting ' + str(len(readList)) + ' .fast5 reads from root directory : ' + self.inputRootDirectory)
          
            self.createSequenceFileOutputs()  

            # Loop through the new reads, extract each one. 
            for readFilename in readList:
                # TODO: Progress bar stuff here.
                
                self.readListOutput.write(readFilename + '\n')
                self.extractSingleMinIONRead(readFilename)

            self.closeOutputFiles()            
          
            # This was a success, reads were extracted.  Lets move on with our lives.
            return True
        
        except Exception:
            # Top Level exception handling like a pro.
            # This is not really doing anything.
            print 'Exception trying to extract a read:'
            print sys.exc_info()
            traceback.print_exc()
            return False

        
        return False
 

    
    def fetchMetaData(self, hdf):
        for metaDataKeyId, metaDataKey in self.metadataLocations.iteritems():        
            try:
                currentAttribute = hdf[metaDataKey].attrs[self.metadataAttributeNames.get(metaDataKeyId)]
                
                print('Fetching Metadata, for this key:' + metaDataKeyId)
                print('Found Attribute:' + currentAttribute)                
                
                previousAttribute = self.metadataAttributes.get(metaDataKeyId)                 
                print('Already Existing Attribute Value:' + str(previousAttribute))                
                
                if(previousAttribute is None):
                    print('Previous attribute is None.  I will assign it.')
                    #self.metadataAttributes.set(metaDataKeyId,currentAttribute)
                    pass
                    #TODO: "set" method isn't working.  I thought dictionaries could do this, maybe something is wrong.
                elif(previousAttribute == currentAttribute):
                    print('The metadata matches with the previous value.  This is good.')
                    pass
                else:
                    print('Metadata mismatch, on attribute ' + metadataKeyId + '.')
                
    
            except KeyError:
                print ('KeyError: metaData Not found: ' + metaDataKeyId + '.  I don\'t care.')
                pass
            except ValueError:
                print ('Value Error.  Investigate.  metaDataKey=' + 
                    metaDataKey + ' , metaDataKeyId=' + metaDataKeyId) 
                traceback.print_exc()
                raise
            except Exception:
                print('Undefined Exception.')
                traceback.print_exc()
                raise
