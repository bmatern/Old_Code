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

# See the file README.MD for how to set up your anaconda environment.

# You should change the following input variables before you run this program:
# The $Rundate variable is included in the name of the output files, 
# You can fill that variable with whatever text you wish to appear there
# Don't use special characters ("*", "/", "\", "$", "?", etc.)

#ReadInputDirectory="/home/minion/MinIONData/2017.03.08.ReadExtractorData/pass"
ReadInputDirectory="/media/ben/D6945F63945F44E1/MinION/Backup_7_596GB/07-01-2017_R9_4_4x_CI_II_LT/downloads/pass"
#ResultsOutputDirectory="/home/minion/MinIONData/2017.03.08.ReadExtractorData/extracts"
ResultsOutputDirectory="/home/minion/MinIONData/Linh/2017.Apr4.DetectDPBCoamplification/Min7000.07012017"
BarcodeFile="/home/minion/MinIONData/2017.03.08.ReadExtractorData/BarcodeSampleMap.txt"
MinimumLength="7000"
MaximumLength="700099999"
Rundate="2017_Jan_1_LinhsReads"

source activate ReadExtractorEnvironment

cd src

#python MinION_Read_Extractor_CL.py -v
#python MinION_Read_Extractor_CL.py --version -i /home/minion/MinIONData/PoretoolsTest/12barcodesR9/pass 
#python MinION_Read_Extractor_CL.py -h
#python MinION_Read_Extractor_CL.py --help

python MinION_Read_Extractor_CL.py \
 --idir=$ReadInputDirectory \
 --odir=$ResultsOutputDirectory \
 --minlen=$MinimumLength \
 --maxlen=$MaximumLength \
 --rundate=$Rundate


source deactivate
