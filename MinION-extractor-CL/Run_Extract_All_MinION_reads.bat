:: This file is part of MinION-extractor-CL.
::
:: MinION-extractor-CL is free software: you can redistribute it and/or modify
:: it under the terms of the GNU Lesser General Public License as published by
:: the Free Software Foundation, either version 3 of the License, or
:: (at your option) any later version.
::
:: MinION-extractor-CL is distributed in the hope that it will be useful,
:: but WITHOUT ANY WARRANTY; without even the implied warranty of
:: MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
:: GNU Lesser General Public License for more details.
::
:: You should have received a copy of the GNU Lesser General Public License
:: along with MinION-extractor-CL. If not, see <http://www.gnu.org/licenses/>.

:: Version 1.0

:: See the file README.MD for how to set up your anaconda environment.

set ReadInputDirectory="C:\MinION-extractor-CL\pass"
set ResultsOutputDirectory="C:\MinION-extractor-CL\extracts"
set BarcodeFile="C:\MinION-extractor-CL\BarcodeSampleMap.txt"
set MinimumLength="1"
set MaximumLength="999999999"
set Rundate="2017_March_14"

activate ReadExtractorEnvironment &&^
cd src &&^
python MinION_Read_Extractor_CL.py ^
 --idir=%ReadInputDirectory% ^
 --odir=%ResultsOutputDirectory% ^
 --barcode=%BarcodeFile% 
 --minlen=%MinimumLength% ^
 --maxlen=%MaximumLength% ^
 --rundate=%Rundate% &&^
deactivate &&^
PAUSE


