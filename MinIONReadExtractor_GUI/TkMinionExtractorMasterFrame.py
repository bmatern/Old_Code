import os
from os import listdir
from os.path import isfile, join, split
import Tkinter, Tkconstants, tkFileDialog, tkMessageBox

#from TkMinionExtractorSubFrame import TkMinionExtractorSubFrame
#from TkSampleIDDialog import TkSampleIDDialog

# This class represents a GUI element for selecting an input folder.
# It will create several instances of 
# TkMinionExtractorSubFrame
# Which are responsible for specifying output files.
class TkMinionExtractorMasterFrame:
    def __init__(self, master):

        master.title("Fast5 Extractor")

        # This is the directory the python executable is running from.
        FileAndPath = os.path.abspath(__file__)
        self.idir, self.ifile = os.path.split(FileAndPath)

        # GUI options
        self.label_opt = {'fill': Tkconstants.BOTH, 'padx': 10, 'pady': 10}
        self.button_opt = {'fill': Tkconstants.BOTH, 'padx': 50, 'pady': 15}
        self.dir_opt = {'initialdir': self.idir, 'mustexist': True, 
            'parent': master, 'title': 'Choose a directory containing .fast5 files'}

        self.master = master
        self.frame = Tkinter.Frame(self.master)

        self.instructionText = Tkinter.StringVar()
        #self.defaultMessage()
        Tkinter.Label(self.frame, width=80, height=5, textvariable=self.instructionText).pack(self.label_opt)
        self.submitButton=Tkinter.Button(self.frame, text='Choose .fast5 Directory:', command=self.getFast5Directory)
        self.submitButton.focus_set()

        # This is what happens when user hits "Enter" on the keyboard.  Simple stuff.
        def returnEvent(event):
            self.getFast5Directory()
        self.master.bind('<Return>', returnEvent)

        self.submitButton.pack(self.button_opt)

        self.frame.pack()


    def defaultMessage(self):
        self.instructionText.set('Choose a directory.  Either:\n\n'
            + '1) A folder containing .fast5 files\n\n'
            + '2) A parent of multiple .fast5 folders (ex. a \"pass\" folder)')

    # getFast5Directory method is called when the user presses the input directory button
    def getFast5Directory(self):
        self.submitButton.config(state='disabled')
        self.setInputDir(tkFileDialog.askdirectory(**self.dir_opt))
        self.submitButton.config(state='normal')

    # setInputDir is a subprocess to be called when an input directory is selected.
    def setInputDir(self, currentFast5Dir):
        currentFast5Dir = os.path.normpath(currentFast5Dir)

        # The default output directory is called 'extracted_reads'.
        # This folder should be at the same directory tree level as the input folder,
        # To allow the barcoded reads to be in the same directory.
        outputDirectory = os.path.normpath(os.path.join(
            os.path.normpath(currentFast5Dir + os.sep + os.pardir)
            ,'extracted_reads'))

        # Open an extractor for the root directory. 
        # This extractor is run whether or not there are subfolders. 
        self.instructionText.set('Extracting sequences from this root directory:\n'
            + currentFast5Dir)         
        tempDir, tempFileName = os.path.split(currentFast5Dir)  
        #sampleID = self.getSampleID(tempFileName) 
        #if not (sampleID == 'SKIPFOLDER'):      
        #    self.openSubFrame(currentFast5Dir, outputDirectory, sampleID)
        self.defaultMessage()
                    
