"""
Tools for making netCDF files desired by WA State DNR for modeling results.
"""

import netCDF4
import time
import os
import numpy

def make_nc_input(fname_nc, fg, force=False, verbose=True):

    
    if os.path.isfile(fname_nc):
        if force and verbose:
            print('Overwriting ', fname_nc)
        elif not force:
            print('*** netCDF file already exists, \n'\
                + '*** NOT overwriting '\
                + '--- use force==True to overwrite' )
            return -1
    
    with netCDF4.Dataset(fname_nc, 'w') as rootgrp:

        rootgrp.description = "fgmax data for " + fg.id
        rootgrp.history = "Created with input data " + time.ctime(time.time())
        rootgrp.history += " in %s;  " % os.getcwd()
            
        if fg.X is not None:
            x = fg.X[0,:]
            lon = rootgrp.createDimension('lon', len(x))
            longitudes = rootgrp.createVariable('lon','f8',('lon',))
            longitudes[:] = x
            longitudes.units = 'degrees_east'
        else:
            if verbose: print('fg.X is None, not adding x')
            
        if fg.Y is not None:
            y = fg.Y[:,0]
            lat = rootgrp.createDimension('lat', len(y))
            latitudes = rootgrp.createVariable('lat','f8',('lat',))
            latitudes[:] = y
            latitudes.units = 'degrees_north'
        else:
            if verbose: print('fg.Y is None, not adding y')
            
        if fg.Z is not None:
            Z = rootgrp.createVariable('Z','f4',('lat','lon',))
            Z[:,:] = fg.Z.data  # include points that are not fg
            Z.units = 'meters'
        else:
            if verbose: print('fg.Z is None, not adding')
            
        if fg.fgmax_point is not None:
            fgmax_point_var = \
                rootgrp.createVariable('fgmax_point','u1',('lat','lon',))
            fgmax_point_var[:,:] = fg.fgmax_point
        else:
            if verbose: print('fg.fgmax_point is None, not adding')
            
        if fg.force_dry_init is not None:
            force_dry_init = \
                rootgrp.createVariable('force_dry_init','u1',('lat','lon',))
            force_dry_init[:,:] = fg.force_dry_init
        else:
            if verbose: print('fg.force_dry_init is None, not adding')  

        print('Created %s' % fname_nc)            
        if verbose:
            print('History:  ', rootgrp.history) 
        return 0     
        
def write_nc_output(fname_nc, fg, new=False, force=False, 
                    outdir='Unknown', verbose=True):

    from clawpack.clawutil.data import ClawData 
    
    fv = -9999.   # fill_value for netcdf4
    
    if new:
        # first create a new .nc file with X,Y,fgmax_point,force_dry_init:
        result = make_nc_input(fname_nc, fg, force=force, verbose=verbose)
        if result == -1:
            print('*** make_nc_input failed, not appending output')
            return        
        
    if outdir is 'Unknown':
        # Cannot determine tfinal or run_finished time
        tfinal = fv
        run_finished = 'Unknown'
    else:
        claw = ClawData()
        claw.read(outdir+'/claw.data', force=True)

        try:
            if claw.output_style==1:
                tfinal = claw.tfinal
            elif claw.output_style==2:
                tfinal = numpy.array(claw.output_times).max()
        except:
            tfinal = fv
        
        try:
            mtime = os.path.getmtime(outdir+'/timing.txt')
            run_finished = time.ctime(mtime) 
        except:
            run_finished = 'Unknown'
            
    # add fgmax output results to existing file
    print(os.getcwd())
    with netCDF4.Dataset(fname_nc, 'a') as rootgrp:
        if verbose:
            print('Appending data from fg to nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('        fg.id: ', fg.id)
        
        h = rootgrp.variables.get('h', None)
        if (h is not None) and (not force):
            print('*** netCDF file already contains output,\n'\
                + '*** NOT overwriting '\
                + '--- use force==True to overwrite' )
            return
                
        x = numpy.array(rootgrp.variables['lon'])
        y = numpy.array(rootgrp.variables['lat'])
        X,Y = numpy.meshgrid(x,y)
        Z = numpy.array(rootgrp.variables['Z'])
        fgmax_point = numpy.array(rootgrp.variables['fgmax_point'])
        bounding_box = [x.min(),x.max(),y.min(),y.max()]
        
        dx = x[1]-x[0]
        Xclose = numpy.allclose(fg.X, X, atol=0.1*dx)
        Yclose = numpy.allclose(fg.Y, Y, atol=0.1*dx)
        
        if (fg.X.shape != X.shape):
            # for now raise an exception, might want to extent to allow
            # filling only part of input arrays
            print('*** Mismatch of fg with data in nc file:')
            print('fg.X.shape = ',fg.X.shape)
            print('nc  X.shape = ',X.shape)
            print('fg.bounding_box = ',fg.bounding_box())
            print('nc  bounding_box = ',bounding_box)
            raise ValueError('*** Mismatch of fg with data in nc file')
    
        Xclose = numpy.allclose(fg.X, X, atol=0.1*dx)
        Yclose = numpy.allclose(fg.Y, Y, atol=0.1*dx)
        if (not (Xclose and Yclose)):
            raise ValueError('*** Mismatch of fg.X or fg.Y with data in nc file')
            

        rootgrp.history += "Added output " + time.ctime(time.time())
        rootgrp.history += " in %s;  " % os.getcwd()
        
        rootgrp.tfinal = tfinal
        rootgrp.outdir = os.path.abspath(outdir)
        rootgrp.run_finished = run_finished
        
        fgmax_point = rootgrp.variables.get('fgmax_point', None)

        if fg.dz is not None:
            try:
                dz = rootgrp.variables['dz']
            except:
                dz = rootgrp.createVariable('dz','f4',('lat','lon',),
                                            fill_value=fv)
            dz[:,:] = fg.dz
            dz.units = 'meters'
            if verbose: print('    Adding fg.dz to nc file')
        else:
            if verbose: print('fg.dz is None, not adding')

        if fg.B is not None:
            try:
                B = rootgrp.variables['B']
            except:
                B = rootgrp.createVariable('B','f4',('lat','lon',),
                                            fill_value=fv)
            B[:,:] = fg.B
            B.units = 'meters'
            if verbose: print('    Adding fg.B to nc file')
        else:
            if verbose: print('fg.B is None, not adding')
                        
        if fg.h is not None:
            try:
                h = rootgrp.variables['h']
            except:
                h = rootgrp.createVariable('h','f4',('lat','lon',),
                                            fill_value=fv)
            h[:,:] = fg.h
            h.units = 'meters'
            if verbose: print('    Adding fg.h to nc file')
        else:
            if verbose: print('fg.h is None, not adding')
            
        if fg.s is not None:        
            try:
                s = rootgrp.variables['s']
            except:
                s = rootgrp.createVariable('s','f4',('lat','lon',),
                                            fill_value=fv)
            s[:,:] = fg.s
            s.units = 'meters/second'
            if verbose: print('    Adding fg.s to nc file')
        else:
            if verbose: print('fg.s is None, not adding')
            
        if fg.hss is not None:        
            try:
                hss = rootgrp.variables['hss']
            except:
                hss = rootgrp.createVariable('hss','f4',('lat','lon',),
                                            fill_value=fv)
            hss[:,:] = fg.hss
            hss.units = 'meters^3/sec^2'
            if verbose: print('    Adding fg.hss to nc file')
        else:
            if verbose: print('fg.hss is None, not adding')
            
        if fg.hmin is not None:        
            try:
                hmin = rootgrp.variables['hmin']
            except:
                hmin = rootgrp.createVariable('hmin','f4',('lat','lon',),
                                            fill_value=fv)
            # negate hmin so that it is minimum flow depth min(h):
            hmin[:,:] = -fg.hmin
            hmin.units = 'meters'
            if verbose: print('    Adding fg.hmin to nc file')
        else:
            if verbose: print('fg.hmin is None, not adding')
            
        if fg.arrival_time is not None:        
            try:
                arrival_time = rootgrp.variables['arrival_time']
            except:
                arrival_time = rootgrp.createVariable('arrival_time','f4',('lat','lon',),
                                            fill_value=fv)
            arrival_time[:,:] = fg.arrival_time
            arrival_time.units = 'seconds'
            if verbose: print('    Adding fg.arrival_time to nc file')
        else:
            if verbose: print('fg.arrival_time is None, not adding')
            
        print('Created %s' % fname_nc)
        if verbose:
            print('History:  ', rootgrp.history)
            print('\nMetadata:')
            print('  outdir:  ', rootgrp.outdir)
            print('  run_finished:  ', rootgrp.run_finished)
            print('  tfinal:  ', rootgrp.tfinal)

def read_nc(fname_nc, verbose=True):

    import netCDF4
    import time
    import os
    from clawpack.geoclaw import fgmax_tools
                
    def get_as_array(var, fgvar=None):
        if fgvar is None:
            fgvar = var
        a = rootgrp.variables.get(var, None)
        if a is not None:
            if verbose: print('    Loaded %s as fg.%s' % (var,fgvar))
            return numpy.array(a)
        else:
            if verbose: print('    Did not find %s for fg.%s' \
                                % (var,fgvar))
            return None
                    
    fg = fgmax_tools.FGmaxGrid()

    with netCDF4.Dataset(fname_nc, 'r') as rootgrp:
        if verbose:
            print('Reading data to fg from nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('History:  ', rootgrp.history)


                
        x = get_as_array('lon','x')
        y = get_as_array('lat','y')
        
        if (x is None) or (y is None):
            print('*** Could not create grid')
            return None
            
        X,Y = numpy.meshgrid(x,y)
        fg.X = X
        fg.Y = Y
        if verbose:
            print('    Constructed fg.X and fg.Y')
        
        fg.Z = get_as_array('Z')
        fg.B = get_as_array('B')
        fg.fgmax_point = get_as_array('fgmax_point') 
        fg.dz = get_as_array('dz')
        fg.h = get_as_array('h')
        fg.s = get_as_array('s')
        fg.hss = get_as_array('hss')
        fg.hmin = get_as_array('hmin')
        fg.arrival_time = get_as_array('arrival_time')
        fg.force_dry_init = get_as_array('force_dry_init')
        
    if verbose:
        print('Returning FGmaxGrid object fg')
    return fg
