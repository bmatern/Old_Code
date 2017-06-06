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
from Bio import SeqIO
from StringIO import StringIO
import os
from os import listdir
from os.path import isfile, join, split
import sys
import subprocess
import MinIONRead
import pylab
from matplotlib import collections

# This script can be used to extract the fastq and fasta sequences 
# from a directory of fast5 files, from MinION Output
# It expects a single parameter, to specify a search directory.  
# Output files are generated in a read_extracts folder under the input directory.




class MinION_Read_Extractor:
    
    # These keys represent the locations of fastq data within a fast5 file.  
    # This might change, depending on Oxford Nanopore file format changes.
    fastqKeys = {
         'TwoDir'           : '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'
        ,'OneDirTemplate'   : '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
        ,'OneDirComplement' : '/Analyses/Basecall_1D_000/BaseCalled_complement/Fastq'
        ,'Barcoding'        : '/Analyses/Barcoding_000/Barcoding/Fastq'
    }
    
    # These are the locations of the metadata attributes on the Fast5 files.
    metadataKeys=MinIONRead.MetaDataDictionary()
    metadataKeys.set('BarcodingChimaeraVersion',     '/Analyses/Barcoding_000')
    metadataKeys.set('BarcodingDragonetVersion',     '/Analyses/Barcoding_000')
    metadataKeys.set('Basecall1DChimaeraVersion',    '/Analyses/Basecall_1D_000')
    metadataKeys.set('Basecall1DDragonetVersion',    '/Analyses/Basecall_1D_000')
    metadataKeys.set('Basecall2DChimaeraVersion',    '/Analyses/Basecall_2D_000')
    metadataKeys.set('Basecall2DDragonetVersion',    '/Analyses/Basecall_2D_000')
    metadataKeys.set('CalibrationChimaeraVersion',   '/Analyses/Calibration_Strand_000')
    metadataKeys.set('CalibrationDragonetVersion',   '/Analyses/Calibration_Strand_000')
    metadataKeys.set('HairpinSplitChimaeraVersion',  '/Analyses/Hairpin_Split_000')
    metadataKeys.set('HairpinSplitDragonetVersion',  '/Analyses/Hairpin_Split_000')
    metadataKeys.set('EventDetectionMinKNOWVersion', '/Analyses/EventDetection_000')
    metadataKeys.set('DeviceID',                     '/UniqueGlobalKey/tracking_id')
    metadataKeys.set('ProtocolsVersionName',         '/UniqueGlobalKey/tracking_id')
    metadataKeys.set('TrackingVersion',              '/UniqueGlobalKey/tracking_id')
    
    
    # These are the names of the attributes of the metadata
    metadataAttributeNames=MinIONRead.MetaDataDictionary()
    metadataAttributeNames.set('BarcodingChimaeraVersion',     'chimaera version')
    metadataAttributeNames.set('BarcodingDragonetVersion',     'dragonet version')
    metadataAttributeNames.set('Basecall1DChimaeraVersion',    'chimaera version')
    metadataAttributeNames.set('Basecall1DDragonetVersion',    'dragonet version')
    metadataAttributeNames.set('Basecall2DChimaeraVersion',    'chimaera version')
    metadataAttributeNames.set('Basecall2DDragonetVersion',    'dragonet version')
    metadataAttributeNames.set('CalibrationChimaeraVersion',   'chimaera version')
    metadataAttributeNames.set('CalibrationDragonetVersion',   'dragonet version')
    metadataAttributeNames.set('HairpinSplitChimaeraVersion',  'chimaera version')
    metadataAttributeNames.set('HairpinSplitDragonetVersion',  'dragonet version')
    metadataAttributeNames.set('EventDetectionMinKNOWVersion', 'version')
    metadataAttributeNames.set('DeviceID',                     'device_id')
    metadataAttributeNames.set('ProtocolsVersionName',         'protocols_version_name')
    metadataAttributeNames.set('TrackingVersion',              'version')
    
    # Display progress text every x reads.  250 is a good choice?
    progressBarBlockSize=250
    
    inputRootDirectory    = None
    outputResultDirectory = None
    minimumReadLength     = 0
    maximumReadLength     = 0
    minimumQuality        = 0
    barcodeSampleMap      = None
    rundate               = None
        
    def __init__(self, inputRootDirectory=None, outputResultDirectory=None):
        self.inputRootDirectory    =inputRootDirectory
        self.outputResultDirectory =outputResultDirectory

    # Print a directory's worth of reads to fasta/fastq.
    #directory, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount)
    def outputFastxResults(self, directory, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount):
    
        print ('OutputFastxResults:' + directory)

        # I'm naming the output file after the subdirectory it's in, relative to the root
        # Usually, after this step, relativePath will contain the barcode text, for example "BC01"
        relativePath = directory.replace(self.inputRootDirectory, '').replace('/','').replace('\\','')[0:]
        if(len(relativePath) < 1):
            relativePath = 'extracted'
        
        if(self.barcodeSampleMap is not None):
            #print('Barcode Sample Map Found')
            if(relativePath in self.barcodeSampleMap):
                print('A mapping from barcode to sample was found. I will rename ' + relativePath + ' to ' + self.barcodeSampleMap[relativePath] + '\n')
                relativePath = self.barcodeSampleMap[relativePath]
            else:
                print('No mapping from barcode to sample was found. I will still name the output file after ' + relativePath + '\n')
            
        else:
            print('No Barcode Sample Map Found.  I will output results to the file named' + relativePath)
            
        if(self.rundate is not None and len(self.rundate) > 0):
            print('A rundate was found, I will add ' + self.rundate + ' to the filename.\n')
            relativePath = relativePath + '_' + self.rundate
            
        # This output directory should exist
        if not os.path.exists(self.outputResultDirectory):
            os.makedirs(self.outputResultDirectory)
        
        # Store output files in a dictionary, for easy looping
        
        #fastAFileOutputs = {
        #     'TwoDir'           : open( join( self.outputResultDirectory , ( relativePath  + '_TwoDirReads.fasta'    )),'w')
        #    ,'OneDirTemplate'   : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirTempReads.fasta')),'w')
        #    ,'OneDirComplement' : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirCompReads.fasta')),'w')
        #    ,'Barcoding'        : open( join( self.outputResultDirectory , ( relativePath  + '_BarcodingReads.fasta' )),'w')
        #}
    
        #fastQFileOutputs = {
        #     'TwoDir'           : open( join( self.outputResultDirectory , ( relativePath  + '_TwoDirReads.fastq'    )),'w')
        #    ,'OneDirTemplate'   : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirTempReads.fastq')),'w')
        #    ,'OneDirComplement' : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirCompReads.fastq')),'w')
        #    ,'Barcoding'        : open( join( self.outputResultDirectory , ( relativePath  + '_BarcodingReads.fastq' )),'w')
        #}
        
    
        DirectoryResultsFilename = join( self.outputResultDirectory , (relativePath  + '_Readstats.txt'))
        
        extractedLengthGraphValues  = []
        extractedQualityGraphValues = []
    
        # Loop through each read in the directory
        for zeroBasedIndex, read in enumerate(directoryReads):
            
            extractedLengthGraphValues.append(len(read.sequence))
            extractedQualityGraphValues.append(read.averageQuality)
    
            readIndex = zeroBasedIndex+1
            #Progress text.  This is printed to console.
            if( (readIndex == 1) or (readIndex == len(directoryReads)) or (readIndex % self.progressBarBlockSize)==0 ):
                percentage = str( (100.00 * readIndex) / len(directoryReads) )[:6] + '%'
                msg = ('Writing reads to file: '
                    + '(' + str(readIndex) + '/' + str(len(directoryReads)) + ') ... ' + percentage)
                print(msg)    
    
            # Choose an output file.  It's determined by the Read Type.
            foundReadType = False
            for id, key in fastAFileOutputs.iteritems():            
                if (read.readType==id):
                    foundReadType = True
                    SeqIO.write([read.sequence], key, 'fasta')
    
            for id, key in fastQFileOutputs.iteritems():            
                if (read.readType==id):
                    SeqIO.write([read.sequence], key, 'fastq')
    
            if (not foundReadType):
                print ('ReadType Not Found.  This could be a problem, investigate this!')
                raise Exception('ReadType Not Found:' + read.readType)
          
        # Close output files.
        for id, key in fastAFileOutputs.iteritems():
            key.close()
    
        for id, key in fastQFileOutputs.iteritems():
            key.close()
            
        # Make a scatterplot of quality vs length of extracted reads
        self.createScatterPlot('Extracted Reads Q score vs Read Length' 
            , extractedLengthGraphValues, extractedQualityGraphValues
            , 'Read Length', 'Average Read Quality (Q score)'
            , 'ExtractedReadsQualityVsLength.png')
    
        rejectLengthFastaOutput = open( join( self.outputResultDirectory , ( relativePath  + '_LengthRejectedReads.fasta'    )),'w')
        rejectLengthFastqOutput = open( join( self.outputResultDirectory , ( relativePath  + '_LengthRejectedReads.fastq'    )),'w')
        
        rejectQualityFastaOutput = open( join( self.outputResultDirectory , ( relativePath  + '_QualityRejectedReads.fasta'    )),'w')
        rejectQualityFastqOutput = open( join( self.outputResultDirectory , ( relativePath  + '_QualityRejectedReads.fastq'    )),'w')
        
        
        lengthRejectLengthGraphValues  = []
        lengthRejectQualityGraphValues = []
        
        qualityRejectLengthGraphValues  = []
        qualityRejectQualityGraphValues = []
        # Print the length rejected reads.
        for zeroBasedIndex, read in enumerate(lengthRejectReads):
    
            SeqIO.write([read.sequence], rejectLengthFastqOutput, 'fastq')
            SeqIO.write([read.sequence], rejectLengthFastaOutput, 'fasta')
            
            lengthRejectLengthGraphValues.append(len(read.sequence))
            lengthRejectQualityGraphValues.append(read.averageQuality)
            
        # Print the quality rejected reads.
        for zeroBasedIndex, read in enumerate(qualityRejectReads):
    
            SeqIO.write([read.sequence], rejectQualityFastqOutput, 'fastq')
            SeqIO.write([read.sequence], rejectQualityFastaOutput, 'fasta')
            
            qualityRejectLengthGraphValues.append(len(read.sequence))
            qualityRejectQualityGraphValues.append(read.averageQuality)
    
            
        rejectLengthFastaOutput.close()
        rejectLengthFastqOutput.close() 
        
        rejectQualityFastaOutput.close()
        rejectQualityFastqOutput.close()   
        
        # Make a scatterplot of quality vs length of extracted reads
        self.createScatterPlot('Length-Rejected Reads Q score vs Read Length' 
            , lengthRejectLengthGraphValues, lengthRejectQualityGraphValues
            , 'Read Length', 'Average Read Quality (Q score)'
            , 'LengthRejectedReadsQualityVsLength.png')
        
        # Make a scatterplot of quality vs length of extracted reads
        self.createScatterPlot('Quality-Rejected Reads Q score vs Read Length' 
            , qualityRejectLengthGraphValues, qualityRejectQualityGraphValues
            , 'Read Length', 'Average Read Quality (Q score)'
            , 'QualityRejectedReadsQualityVsLength.png')
        
    
        # Write Metadata and Results
        self.writeResultsOutput(DirectoryResultsFilename, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount)
        
        
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
           
        
        
    # writeResultsOutput method calculates simple stats on the reads, and writes relevant metadata to a file.  
    def writeResultsOutput(self, resultsFileName, reads, lengthRejectReads, qualityRejectReads, totalFast5FileCount):
    
        print('Creating Results File:' + resultsFileName)
        resultsOutput = open(resultsFileName, 'w')
     
        resultsOutput.write('Read Extraction Overall Summary\n\n')
        
        resultsOutput.write('Read Directory:\n' + self.inputRootDirectory + '\n')
        resultsOutput.write(str(self.minimumReadLength) + ' Minimum Read Length\n')
        resultsOutput.write(str(self.maximumReadLength) + ' Maximum Read Length\n')
        resultsOutput.write(str(self.minimumQuality) + ' Minimum Quality\n\n')
        
        totalAnalyzedReadCount = (len(reads) + len(lengthRejectReads) + len(qualityRejectReads))
        
        resultsOutput.write(str(totalFast5FileCount) + ' Fast5 files processed, which contained:\n')
        resultsOutput.write(str(totalAnalyzedReadCount) + ' Reads total\n')    
        resultsOutput.write(str(len(reads)) + ' Reads extracted\n')
        resultsOutput.write(str(len(lengthRejectReads)) + ' Reads rejected for length\n')
        resultsOutput.write(str(len(qualityRejectReads)) + ' Reads rejected for quality\n\n')
        
    
        TwoDirReads = []
        OneDirTempReads = []
        OneDirCompReads = []
        BarcodingReads = []
    
        currentMetadata = MinIONRead.MetaDataDictionary()
        
        
    
        for zeroBasedIndex, read in enumerate(reads):
            
            #Progress text.  This is printed to console.
            readIndex = zeroBasedIndex+1
            if( (readIndex == 1) or (readIndex == len(reads)) or (readIndex % self.progressBarBlockSize)==0 ):
                percentage = str( (100.00 * readIndex) / len(reads) )[:6] + '%'
                msg = ('Compiling Results: '
                    + '(' + str(readIndex) + '/' + str(len(reads)) + ') ... ' + percentage)
                print(msg)    
    
            # Get the read type.
            if (read.readType=='TwoDir'):
                TwoDirReads.append(read)
            elif (read.readType=='OneDirTemplate'):
                OneDirTempReads.append(read)
            elif (read.readType=='OneDirComplement'):
                OneDirCompReads.append(read)
            elif (read.readType=='Barcoding'):
                BarcodingReads.append(read)
            else:
                print ('ReadType Not Found.  This could be a problem...')
    
            # Record Metadata from each read.
            try:
                for id, value in currentMetadata.iteritems():
                    # if we haven't found this metadata yet
                    if (value==''):
                        currentMetadata.set(id , read.get(id))
                    # if the current read matches the value we already have.
                    elif (value==read.get(id)):
                        # I expect this.  Do nothing.
                        pass
                    # if the current read is blank, but we already had a version
                    elif (read.get(id) == ''):
                        # Nothing to do here either, there is no conflict.
                        pass                
                    # apparently the metadata on this read doesn't match what we expect.
                    else:
                        #TODO: This exception was raised once, between TrackingVersion 1.1.20 and 1.1.21.  
                        # I don't know why the format changed, but I cant run it when this exception is raised.
                        # Changing it to a print statement instead....
                        
                        print('Metadata versions don\'t match:' + id + ': ' + value + ' and ' + read.get(id))
                        #raise Exception('Metadata versions don\'t match:' + id + ': ' + value + ' and ' + read.get(id))
     
            except Exception:
                print('Undefined Exception.')
                raise
            
        self.writeReadStatsSub('TwoDir', TwoDirReads, lengthRejectReads, qualityRejectReads,  resultsOutput)
        self.writeReadStatsSub('OneDirTemplate', OneDirTempReads,  lengthRejectReads, qualityRejectReads, resultsOutput)
        self.writeReadStatsSub('OneDirComplement', OneDirCompReads, lengthRejectReads, qualityRejectReads,  resultsOutput)
        self.writeReadStatsSub('Barcoding', BarcodingReads, lengthRejectReads, qualityRejectReads,  resultsOutput)
    
        resultsOutput.write('Metadata:\n')
        for id, value in currentMetadata.iteritems():
            resultsOutput.write(id + ':' + value + '\n')
        resultsOutput.write('\n')
    
        resultsOutput.close()
        
    
        
    # Short method to print out summary stats for a group of reads
    def writeReadStatsSub(self, readType, reads, lengthRejectReads, qualityRejectReads, outputFile):
        
        # Count the rejected reads quick
        lengthRejectCount = 0
        qualityRejectCount = 0
        
        for zeroBasedIndex, read in enumerate(lengthRejectReads):
            if(read.readType == readType):
                lengthRejectCount += 1
                
        for zeroBasedIndex, read in enumerate(qualityRejectReads):
            if(read.readType == readType):
                qualityRejectCount += 1
        
        
        if(len(reads) > 0):
            reads.sort(key=lambda x: len(x.sequence), reverse=False)
            outputFile.write(readType + ' Reads Summary\n' + str(len(reads)) + ' Reads extracted'
            + '\n' + str(lengthRejectCount) + ' Reads rejected for length'
            + '\n' + str(qualityRejectCount) + ' Reads rejected for quality'
            + '\n' + str(len(reads[0].sequence)) + ' Minimum Length'       
            + '\n' + str(len(reads[len(reads) - 1].sequence)) + ' Maximum Length'
            + '\n' + str(sum(len(read.sequence) for read in reads) /len(reads)) + ' Mean Length'
            + '\n' + str(sum(read.averageQuality for read in reads) /len(reads)) + ' Mean Quality'
            + '\n\n')
        
        else:
            outputFile.write(readType + ' Reads Summary\n0 Reads extracted'
            + '\n' + str(lengthRejectCount) + ' Reads rejected for length'
            + '\n' + str(qualityRejectCount) + ' Reads rejected for quality'
            + '\n\n')
    
    
    # Recursive method to extract data from all fast5 files in a directory
    def extractAllFast5(self, directory):
        print('Opening directory ' + directory)
        
        # Get children files and folders    
        fileNames = [f for f in listdir(directory) if isfile(join(directory, f))]
        fileNames.sort()
        subdirectoryNames = [f for f in listdir(directory) if not isfile(join(directory, f))]
        subdirectoryNames.sort()
    
        print(str(len(subdirectoryNames)) + ' subdirectories and ' + str(len(fileNames)) + ' files found.')
    
        #Recurse into each subdirectory.  
        if (len(subdirectoryNames) > 0):
            for subdir in subdirectoryNames:
                if(subdir != 'read_extracts'):
                    self.extractAllFast5(os.path.join(directory, subdir))
    
        #For each file, check all possible fastqKeys for fasta results.  Print them based on the fasta key.
        if (len(fileNames) > 0):
    
            directoryReads = []
            qualityRejectReads = []
            lengthRejectReads = []
            
            totalFast5FileCount = 0
            
            #fast5GoodFileCounter = 0
            #fast5QualityRejectCounter = 0
            #fast5LengthRejectCounter = 0
            
    
            for zeroBasedIndex, fileName in enumerate(fileNames):
    
                readIndex = zeroBasedIndex+1
                #Progress text.  This is printed to console.
                if( (readIndex == 1) or (readIndex == len(fileNames)) or (readIndex % self.progressBarBlockSize)==0 ):
                    percentage = str( (100.00 * readIndex) / len(fileNames) )[:6] + '%'
                    msg = ('Scanning fast5 file: '
                        + '(' + str(readIndex) + '/' + str(len(fileNames)) + ') ... ' + percentage)
                    print(msg)    
    
                combinedFileName = join(directory, fileName)
    
              
                # Only want fast5 files.
                if('fast5' in fileName):    
                    
                    totalFast5FileCount += 1 
                       
                    try:
    
                        hdf = h5py.File(combinedFileName, 'r')
                        for id, key in self.fastqKeys.iteritems():
    
                            #print('id:' + id)
    
                            sequenceFound = False
                            rec = None
                            # Attempt to read the fastq file at this location
                            try:
                                fq = hdf[key][()]
    
                                rec = SeqIO.read(StringIO(fq), "fastq")
                                rec.id += "_" + id + " " + SoftwareVersion
    
                                # Add the location of the read file in the header.  
                                # nanopolish expects to find the file location in the last token of the fasta header.
                                #https://github.com/jts/nanopolish/blob/45335166365ccc46a3dd434f5b937f1594e103f8/src/common/nanopolish_fast5_map.cpp
                                rec.description += " " + combinedFileName
                                sequenceFound = True
    
                            except KeyError:
                                #print('Probably Okay, Exception reading the fastq, key =' + key)
                                # pass is just a dummy statement because these exceptions are completely expected and okay.  
                                # There is no data in this key, and that's fine most of the time.
                                pass
    
                            if (sequenceFound):
                                try:
                                    # "qualities" will contain letter-by-letter quality scores 
                                    # for each sequenced nucleotides. 
                                    # An algebraic mean doesn't makes much sense for logarithmic phred scores.
                                    # Oh well.
                                    qualities = rec.letter_annotations["phred_quality"]
                                    currentAvgQuality = numpy.mean(qualities)
    
                                    currentMinIONRead = MinIONRead.MinIONRead()
                                    currentMinIONRead.id = rec.id
                                    currentMinIONRead.readType = id
                                    currentMinIONRead.sequence = rec
                                    currentMinIONRead.averageQuality = currentAvgQuality
    
                                    # Fetch Metadata.
                                    self.fetchMetaData(hdf, currentMinIONRead)
                                      
                                    # filter read lengths.  Ignore filter if cutoff is 0, the default
                                    # If (len>min or min=0) AND (len<max or max=0)
                                    #if( (len(currentMinIONRead.sequence) >= minimumReadLength or minimumReadLength == 0)
                                    #    and (len(currentMinIONRead.sequence) <= maximumReadLength or maximumReadLength == 0) 
                                    #    ):          
                                    #    directoryReads.append(currentMinIONRead)
                                        
                                    # If the read length is outside the min/max
                                    if( (len(currentMinIONRead.sequence) < self.minimumReadLength and self.minimumReadLength != 0)
                                        or (len(currentMinIONRead.sequence) > self.maximumReadLength and self.maximumReadLength != 0) 
                                        ):          
                                        lengthRejectReads.append(currentMinIONRead)
                                        
                                    # If read is below quality
                                    elif( currentMinIONRead.averageQuality < self.minimumQuality and self.minimumQuality != 0
                                        ):          
                                        qualityRejectReads.append(currentMinIONRead)
                                        
                                    # Else we have a good quality read that we want to extract.
                                    else:
                                        directoryReads.append(currentMinIONRead)
    
                                except Exception, e:
                                    print('Exception occured when reading a .fast5 file: ' + str(e))
                                    raise
                     
                        hdf.close()
                    except Exception, e:
                        print('Exception when opening .fast5 file: ' + fileName + '\nProbably ok if this doesn\'t happen often: ' + str(e))
     
            #Here is where I should print resulting datafiles
            self.outputFastxResults(directory, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount)
            print('Done extracting fast5 from this directory:' + directory)
    
    def fetchMetaData(self, hdf, currentMinIONRead):
        for metaDataKeyId, metaDataKey in self.metadataKeys.iteritems():        
            try:
                currentAttribute = hdf[metaDataKey].attrs[self.metadataAttributeNames.get(metaDataKeyId)]
                currentMinIONRead.set(metaDataKeyId , currentAttribute)
    
            except KeyError:
                #print ('KeyError: metaData Not found: ' + metaDataKeyId + '.  I don\'t care.')
                pass
            except ValueError:
                print ('Value Error.  Investigate.  metaDataKey=' + 
                    metaDataKey + ' , metaDataKeyId=' + metaDataKeyId) 
                raise
            except Exception:
                print('Undefined Exception.')
                raise
