#!/usr/bin/env python

"""
To add to: amrclaw.data

Base AMRClaw data class for writing out data parameter files.

flagregions.data

flagregion attributes:
    minlevel
    maxlevel
    t1
    t2
    spatial_region_type = 1 for rectangle, 2 for RuledRectangle
    if spatial_region_type == 1:
        spatial_region = [x1,x2,y1,y2] 
    if spatial_region_type == 2:
        spatial_region_file = file describing RuledRectangle
    
Possible future attributes:
    flag_method - (e.g. adjoint, wave_tolerance) to use in this region
    tolerance - for flagging in this region
    
"""

from __future__ import print_function
from __future__ import absolute_import
import os

import clawpack.clawutil.data
import six
from six.moves import range
import numpy as np

import region_tools

class FlagRegion(clawpack.clawutil.data.ClawData):
    
    def __init__(self, num_dim, region=None):
        r"""
        Note: for backward compatibility, can pass in a old-style region 
        in the form (in 2d):
            [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
        and it will be converted to the new form.
        
        """
        
        super(FlagRegion,self).__init__()
        self.add_attribute('name','')
        self.add_attribute('num_dim',num_dim)
        self.add_attribute('minlevel',None)
        self.add_attribute('maxlevel',None)
        self.add_attribute('t1',None)
        self.add_attribute('t2',None)
        self.add_attribute('spatial_region_type',None)
        self.add_attribute('spatial_region',None)
        self.add_attribute('spatial_region_file',None)
        
        if region is not None:
            self.convert_old_region(region)

        
    def convert_old_region(region):
        """
        Take a list region = [minlevel, maxlevel, t1, t2, x1, x2, y1, y2]
        in the old style and convert to a new flagregion.
        """
        self.minlevel = region[0]
        self.maxlevel = region[1]
        self.t1 = region[2]
        self.t2 = region[3]
        self.spatial_region_type = 1  # rectangle
        self.spatial_region = region[4:]
        
    def read_spatial_region(self, fname=None):
        """
        Reads a ruled rectangle and assigns to self.spatial_region
        But this isn't really supported yet, since cannot write to 
        flagregions.data.
        """
        if fname is None:
            fname = self.spatial_region_file
        if not os.path.isfile(fname):
            print('*** Could not find file ',fname)
        else:
            self.spatial_region = region_tools.RuledRectangle(fname)
            
        
# ==============================================================================
#  Region data object
class FlagRegionData(clawpack.clawutil.data.ClawData):
    r"""To replace RegionData, allowing more flexibility."""

    def __init__(self,flagregions=None,num_dim=2):

        super(FlagRegionData,self).__init__()

        if flagregions is None or not isinstance(regions,list):
            self.add_attribute('flagregions',[])
        else:
            self.add_attribute('flagregions',regions)
        self.add_attribute('num_dim',num_dim)


    def write(self,out_file='flagregions.data',data_source='setrun.py'):
        """
        Write out the flagregion.data file that will be read by the Fortran
        code.  
        """

        self.open_data_file(out_file,data_source)

        self.data_write(value=len(self.flagregions),alt_name='num_flagregions')
        for flagregion in self.flagregions:
            flagregion._out_file = self._out_file
            # write minlevel,maxlevel as integers, t1,t2 as floats:
            flagregion.data_write()
            flagregion.data_write('name')
            flagregion.data_write('minlevel')
            flagregion.data_write('maxlevel')
            flagregion.data_write('t1')
            flagregion.data_write('t2')
            flagregion.data_write('spatial_region_type')
            
            if flagregion.spatial_region_type == 1:
                # spatial_region is an extent [x1,x2,y1,y2]
                assert len(flagregion.spatial_region) == 4, \
                    "*** FlagRegionData.write requires len(spatial_region) = 4"
                flagregion.data_write('spatial_region')
                
            elif flagregion.spatial_region_type == 2:
                if type(flagregion.spatial_region_file) is str:
                    # assumed to be path to RuledRectangle file
                    flagregion.data_write('spatial_region_file')
                else:
                    raise ValueError('*** FlagRegionData.write requires path to RuledRectangle')
                    
        self.close_data_file()


    def read(self, path='flagregions.data', force=False):
        r"""Read in flagregion data from file at path
        """

        with open(os.path.abspath(path),'r') as data_file:
            # Read past comments and blank lines
            header_lines = 0
            ignore_lines = True
            while ignore_lines:
                line = data_file.readline()
                if line[0] == "#" or len(line.strip()) == 0:
                    header_lines += 1
                else:
                    break
    
            # Read in number of flagregions        
            num_flagregions, tail = line.split("=:")
            num_flagregions = int(num_flagregions)
            varname = tail.split()[0]
            if varname != "num_flagregions":
                raise IOError("It does not appear that this file contains flagregion data.")
    
            print('Will read %i flagregions' % num_flagregions)
            
            # Read in each flagregion
            self.flagregions = []
            for n in range(num_flagregions):
                flagregion = FlagRegion(num_dim=2)
                line = data_file.readline() # blank
                line = data_file.readline().split()
                flagregion.name = str(line[0])[1:-1] # strip extra quotes
                line = data_file.readline().split()
                flagregion.minlevel = int(line[0])
                line = data_file.readline().split()
                flagregion.maxlevel = int(line[0])
                line = data_file.readline().split()
                flagregion.t1 = float(line[0])
                line = data_file.readline().split()
                flagregion.t2 = float(line[0])
                line = data_file.readline().split()
                flagregion.spatial_region_type = int(line[0])
                
                line = data_file.readline().split()
                if flagregion.spatial_region_type == 1:
                    flagregion.spatial_region = [float(val) for val in line[:4]]
                elif flagregion.spatial_region_type == 2:
                    # the next line is assumed to be path to a
                    # RuledRectangle file, strip quotes
                    flagregion.spatial_region_file = line[0][1:-1] 
                    # read the file: 
                    flagregion.read_spatial_region()
                else:
                    raise ValueError('*** Unrecognized spatial_region_type')
                
                self.flagregions.append(flagregion)

