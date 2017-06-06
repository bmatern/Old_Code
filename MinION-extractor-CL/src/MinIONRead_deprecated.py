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

# This emulates a dictionary containing the metadata values.
class MetaDataDictionary:  
    def __init__(self):
        self.metadata={
             'BarcodingChimaeraVersion'     : ''
            ,'BarcodingDragonetVersion'     : ''
            ,'Basecall1DChimaeraVersion'    : ''
            ,'Basecall1DDragonetVersion'    : ''
            ,'Basecall2DChimaeraVersion'    : ''
            ,'Basecall2DDragonetVersion'    : ''
            ,'CalibrationChimaeraVersion'   : ''
            ,'CalibrationDragonetVersion'   : ''
            ,'EventDetectionMinKNOWVersion' : ''
            ,'HairpinSplitChimaeraVersion'  : ''
            ,'HairpinSplitDragonetVersion'  : ''
            ,'DeviceID'                     : ''
            ,'ProtocolsVersionName'         : ''
            ,'TrackingVersion'              : ''
        }

    # These methods make my life easier and make this class look like a dictionary.
    # Perhaps you are smarter than I am.
    # TODO: This is silly, not sure why I'm storing metadata like this. Cant it just be a class variable?
    def get(self, index):
        return self.metadata.get(index)

    def set(self, index, value):
        self.metadata[index] = value

    def iteritems(self):
        return self.metadata.iteritems()

# Information about a MinION Read
class MinIONRead:
    def __init__(self):
        self.id=None
        self.readType=None
        self.barcode=None
        self.sequence=None
        self.averageQuality=0
        self.filePath=None
        self.metadata=MetaDataDictionary()

    def get(self, index):
        return self.metadata.get(index)

    def set(self, index, value):
        self.metadata.set(index,value)


