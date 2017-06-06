#MinION-extractor-CL
A Commandline tool to extract reads in .fasta and .fastq format, and relevant metadata from the .fast5 files created by a MinION Sequencing run.

##Download the executable
When the compiled executable is available, i'll put a link for the file download here.

##Run using Python
Run this program using standard python call:  
```
python MinION_Read_Extractor_CL.py /home/minion/MinIONData/SampleData
```
The script expects a single parameter, to specify the input directory.

##To configure Anaconda
Install Anaconda for python 2.7.  
https://www.continuum.io/downloads  
To set up the environment in anaconda:  

Linux/Mac:  
```
conda create --name ReadExtractorEnvironment biopython h5py six  
source activate ReadExtractorEnvironment  
pip install pyinstaller packaging  
source deactivate  
```  
Windows:  
```  
conda create --name ReadExtractorEnvironment biopython h5py six pywin32  
call activate ReadExtractorEnvironment && pip install pyinstaller packaging && call deactivate  
```

##Run using a bash or .bat script using Anaconda
You can execute the following scripts to run the tool inside of Anaconda:  
Linux/Mac:  
```
bash Run_Extract_All_MinION_reads.sh  
```  
Windows:  
```  
Run_Extract_All_MinION_reads.bat  
```  

##Input Data
This tool has no user interface.  The user will specify the input folder in a single parameter.  
You may specify a folder containing fast5 reads, or alternatively, you can specify a parent folder containing any number of subfolders containing barcoded or sorted samples. 

##Output Data
Reads are extracted in .fasta and .fastq format.  When the input folder is selected, the MinION-extractor-GUI will specify a default output directory directly underneath the input directory, called "read_extracts"  This tool will recursively parse subdirectories and output data files named by directory structure.  
















