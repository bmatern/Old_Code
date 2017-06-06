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

import os
from os.path import join
from Bio import SeqIO

# Information about a MinION Read Output file
# This class should act like a file
# I should learn inheritance the proper way.
class MinIONReadOutputFile:
    def __init__(self, rundate, barcode, readType, outputDirectory, fileFormat):
   
        self.barcode = barcode
        self.readType = readType
        self.fileFormat = fileFormat
        
        shortFileName = str(rundate) + '_' + str(barcode) + '_' + str(readType) + '_reads.' + str(fileFormat) 
        self.fileName = os.path.join(outputDirectory, shortFileName)
        
        # Just open up the file for appending, I'll close it later.
        self.outputFile = open(self.fileName ,'a')        
          
        
    def write(self, text):
        if self.outputFile.closed:
            self.open()
        self.outputFile.write(text)
    
    def close(self):
        if not self.outputFile.closed:
            self.outputFile.close()
        
    def open(self):
        if self.outputFile.closed:
            self.outputFile = open(self.fileName ,'a')  
            
    def getFastqReads(self):
        self.close()
        fastqReads = SeqIO.parse(self.fileName, 'fastq')
        return fastqReads
            
       
