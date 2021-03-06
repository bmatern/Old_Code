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

SoftwareVersion = "MinION-extractor-GUI Version 1.0"

import Tkinter
import sys
from TkMinionExtractorMasterFrame import TkMinionExtractorMasterFrame

# This is the main entry point for the 2d Read extractor GUI.
# This script can be used to extract the fastq or fasta read sequences 
# from a directory of fast5 files, generated by MinION           
if __name__=='__main__':

        # No parameters are expected at all.  sys.argv[0] doesn't count.
        if (len(sys.argv) == 1):
            print('\n***Running the MinION 2D Read Extractor GUI***\n')

            root = Tkinter.Tk()
            app = TkMinionExtractorMasterFrame(root)
            root.mainloop()

            print('Done.  Yay.')

        # Print the Software Version
        #elif (len(sys.argv) == 2 and (
        #    sys.argv[1].lower() == '-v' or 
        #    sys.argv[1].lower() == '--version' or 
        #    sys.argv[1].lower() == '-version')        
        #):
        #    print (SoftwareVersion)

        # You called the software wrong.  
        #else:
        #    print("usage:\n" + 
        #        "\tRun this program using standard python call:\n" + 
        #        "\t$python Extract2dReadsMain.py\n" + 
        #        "\tbiopython and hdf5 must be accessible in your python environment.\n"
        #        "\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
        #    )

