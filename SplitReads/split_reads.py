import sys
import os
import random
import math
from Bio import SeqIO
from Bio import SeqRecord

def splitReads():

    inputFilename = sys.argv[1]
    tokens = os.path.splitext(inputFilename)
    outputFilename = tokens[0] + '.split' + tokens[1]
    print ('Input file:' + inputFilename)
    print ('Output file:' + outputFilename)

    #print ('arg0:' + (os.path.splitext(inputFilename)[0]))
    #print ('arg1:' + (os.path.splitext(inputFilename)[1]))
    #print ('arg2:' + (os.path.splitext(inputFilename)[2]))
    records = SeqIO.parse(inputFilename, "fasta")

    resultsOutput = open(outputFilename, 'w')

    minimumLength = 2000

    #print (str(len(records)) + ' records found.')
    print ('Splitting reads...')

    for index, record in enumerate(records):
        #print(str(index) + ':' + str(record.id) + ' ' + str(len(record)) + ' ' + str(len(record.seq)))

        if(len(record.seq) > minimumLength * 3):

            splitIndex = int(math.floor(random.random() * len(record)))
            if (splitIndex < minimumLength):
                splitIndex += minimumLength

            if(splitIndex > (len(record) - minimumLength)):
                splitIndex -= minimumLength
                pass

            sequence1 = record.seq[0:splitIndex]
            #print ('SEQUENCE1:' + sequence1 )
            sequence2 = record.seq[splitIndex:len(record.seq)]

            #sequence1Rec = SeqIO.parse(sequence1, "fasta").next()
            #sequence2Rec = SeqIO.parse(sequence2, "fasta").next()

            sequence1Rec = record.upper()
            sequence2Rec = record.upper()

            sequence1Rec.seq = sequence1
            sequence2Rec.seq = sequence2

            sequence1Rec.id = record.id + '_SplitRead1'
            sequence2Rec.id = record.id + '_SplitRead2'

            SeqIO.write([sequence1Rec], resultsOutput, 'fasta')
            SeqIO.write([sequence2Rec], resultsOutput, 'fasta')


            #print (str(len(sequence1)) + ' + ' + str(len(sequence2)) + ' = ' + str(len(record.seq)))

            
            #print('features:' + str(len(record.features)))
            #print ('random:' + str(random.random()) + '...' + str(splitIndex))
            
            pass

        elif(len(record.seq) > 2000):
            # 2000 < length <= 6000
            #add this read by itself
            SeqIO.write([record], resultsOutput, 'fasta')
            pass
        else:
            #Discard this record, do nothing.  It's toast.
            pass

    resultsOutput.close()

if __name__=='__main__':
    try:
        splitReads()
        print('Done.  Yay.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise
