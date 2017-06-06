    # Recursive method to extract data from all fast5 files in a rootDirectory
    def extractAllFast5_old_delete(self, directory):
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
                                    print sys.exc_info()
                                    traceback.print_exc()
                                    raise
                     
                        hdf.close()
                    except Exception, e:
                        print('Exception when opening .fast5 file: ' + fileName + '\nProbably ok if this doesn\'t happen often:\n' + str(e))
                        print sys.exc_info()
                        traceback.print_exc()
     
            #Here is where I should print resulting datafiles
            self.outputFastxResults(directory, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount)
            print('Done extracting fast5 from this directory:' + directory)
            
            
            
            
            
            
    # Print a directory's worth of reads to fasta/fastq.
    #directory, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount)
    def outputFastxResults_old_delete(self, directory, directoryReads, lengthRejectReads, qualityRejectReads, totalFast5FileCount):
    
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
        fastAFileOutputs = {
             'TwoDir'           : open( join( self.outputResultDirectory , ( relativePath  + '_TwoDirReads.fasta'    )),'w')
            ,'OneDirTemplate'   : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirTempReads.fasta')),'w')
            ,'OneDirComplement' : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirCompReads.fasta')),'w')
            ,'Barcoding'        : open( join( self.outputResultDirectory , ( relativePath  + '_BarcodingReads.fasta' )),'w')
        }
    
        fastQFileOutputs = {
             'TwoDir'           : open( join( self.outputResultDirectory , ( relativePath  + '_TwoDirReads.fastq'    )),'w')
            ,'OneDirTemplate'   : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirTempReads.fastq')),'w')
            ,'OneDirComplement' : open( join( self.outputResultDirectory , ( relativePath  + '_OneDirCompReads.fastq')),'w')
            ,'Barcoding'        : open( join( self.outputResultDirectory , ( relativePath  + '_BarcodingReads.fastq' )),'w')
        }
        
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
      