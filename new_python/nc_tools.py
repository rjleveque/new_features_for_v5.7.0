
            

def make_nc_input(fname_nc, fgm, force=False, verbose=True):

    import netCDF4
    import time
    import os
    
    if os.path.isfile(fname_nc):
        if force and verbose:
            print('Overwriting ', fname_nc)
        elif not force:
            print('*** netCDF file already exists, \n'\
                + '*** NOT overwriting '\
                + '--- use force==True to overwrite' )
            return -1
    
    with netCDF4.Dataset(fname_nc, 'w') as rootgrp:

        rootgrp.description = "fgmax data for " + fgm.id
        rootgrp.history = "Created with input data " + time.ctime(time.time())
        rootgrp.history += " in %s;  " % os.getcwd()
            
        if fgm.X is not None:
            x = fgm.X[0,:]
            lon = rootgrp.createDimension('lon', len(x))
            longitudes = rootgrp.createVariable('lon','f8',('lon',))
            longitudes[:] = x
            longitudes.units = 'degrees_east'
        else:
            if verbose: print('fgm.X is None, not adding x')
            
        if fgm.Y is not None:
            y = fgm.Y[:,0]
            lat = rootgrp.createDimension('lat', len(y))
            latitudes = rootgrp.createVariable('lat','f8',('lat',))
            latitudes[:] = y
            latitudes.units = 'degrees_north'
        else:
            if verbose: print('fgm.Y is None, not adding y')
            
        if fgm.Z is not None:
            Z = rootgrp.createVariable('Z','f4',('lat','lon',))
            Z[:,:] = fgm.Z.data  # include points that are not fgmax_points
            Z.units = 'meters'
        else:
            if verbose: print('fgm.Z is None, not adding')
            
        if fgm.fgmax_point is not None:
            fgmax_point_var = \
                rootgrp.createVariable('fgmax_point','u1',('lat','lon',))
            fgmax_point_var[:,:] = fgm.fgmax_point
        else:
            if verbose: print('fgm.fgmax_point is None, not adding')
            
        if fgm.force_dry_init is not None:
            force_dry_init = \
                rootgrp.createVariable('force_dry_init','u1',('lat','lon',))
            force_dry_init[:,:] = fgm.force_dry_init
        else:
            if verbose: print('fgm.force_dry_init is None, not adding')  

        print('Created %s' % fname_nc)            
        if verbose:
            print('History:  ', rootgrp.history) 
        return 0     
        
def write_nc_output(fname_nc, fgm, new=False, force=False, 
                    outdir='Unknown', verbose=True):

    import netCDF4
    import time
    import os
    from clawpack.clawutil.data import ClawData 
    
    fv = -9999.   # fill_value for netcdf4
    
    if new:
        # first create a new .nc file with X,Y,fgmax_point,force_dry_init:
        result = make_nc_input(fname_nc, fgm, force=force, verbose=verbose)
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
        tfinal = claw.tfinal
        
        try:
            mtime = os.path.getmtime(outdir+'/timing.txt')
            run_finished = time.ctime(mtime) 
        except:
            run_finished = 'Unknown'
            
    # add fgmax output results to existing file
    with netCDF4.Dataset(fname_nc, 'a') as rootgrp:
        if verbose:
            print('Appending data from fgm to nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('        fgm.id: ', fgm.id)
        
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
        Xclose = numpy.allclose(fgm.X, X, atol=0.1*dx)
        Yclose = numpy.allclose(fgm.Y, Y, atol=0.1*dx)
        
        if (fgm.X.shape != X.shape):
            # for now raise an exception, might want to extent to allow
            # filling only part of input arrays
            print('*** Mismatch of fgm with data in nc file:')
            print('fgm.X.shape = ',fgm.X.shape)
            print('nc  X.shape = ',X.shape)
            print('fgm.bounding_box = ',fgm.bounding_box())
            print('nc  bounding_box = ',bounding_box)
            raise ValueError('*** Mismatch of fgm with data in nc file')
    
        Xclose = numpy.allclose(fgm.X, X, atol=0.1*dx)
        Yclose = numpy.allclose(fgm.Y, Y, atol=0.1*dx)
        if (not (Xclose and Yclose)):
            raise ValueError('*** Mismatch of fgm.X or fgm.Y with data in nc file')
            

        rootgrp.history += "Added output " + time.ctime(time.time())
        rootgrp.history += " in %s;  " % os.getcwd()
        
        rootgrp.tfinal = tfinal
        rootgrp.outdir = os.path.abspath(outdir)
        rootgrp.run_finished = run_finished
        
        fgmax_point = rootgrp.variables.get('fgmax_point', None)

        if fgm.dz is not None:
            try:
                dz = rootgrp.variables['dz']
            except:
                dz = rootgrp.createVariable('dz','f4',('lat','lon',),
                                            fill_value=fv)
            dz[:,:] = fgm.dz
            dz.units = 'meters'
            if verbose: print('    Adding fgm.dz to nc file')
        else:
            if verbose: print('fgm.dz is None, not adding')

        if fgm.B is not None:
            try:
                B = rootgrp.variables['B']
            except:
                B = rootgrp.createVariable('B','f4',('lat','lon',),
                                            fill_value=fv)
            B[:,:] = fgm.B
            B.units = 'meters'
            if verbose: print('    Adding fgm.B to nc file')
        else:
            if verbose: print('fgm.B is None, not adding')
                        
        if fgm.h is not None:
            try:
                h = rootgrp.variables['h']
            except:
                h = rootgrp.createVariable('h','f4',('lat','lon',),
                                            fill_value=fv)
            h[:,:] = fgm.h
            h.units = 'meters'
            if verbose: print('    Adding fgm.h to nc file')
        else:
            if verbose: print('fgm.h is None, not adding')
            
        if fgm.s is not None:        
            try:
                s = rootgrp.variables['s']
            except:
                s = rootgrp.createVariable('s','f4',('lat','lon',),
                                            fill_value=fv)
            s[:,:] = fgm.s
            s.units = 'meters/second'
            if verbose: print('    Adding fgm.s to nc file')
        else:
            if verbose: print('fgm.s is None, not adding')
            
        if fgm.hss is not None:        
            try:
                hss = rootgrp.variables['hss']
            except:
                hss = rootgrp.createVariable('hss','f4',('lat','lon',),
                                            fill_value=fv)
            hss[:,:] = fgm.hss
            hss.units = 'meters^3/sec^2'
            if verbose: print('    Adding fgm.hss to nc file')
        else:
            if verbose: print('fgm.hss is None, not adding')
            
        if fgm.hmin is not None:        
            try:
                hmin = rootgrp.variables['hmin']
            except:
                hmin = rootgrp.createVariable('hmin','f4',('lat','lon',),
                                            fill_value=fv)
            # negate hmin so that it is minimum flow depth min(h):
            hmin[:,:] = -fgm.hmin
            hmin.units = 'meters'
            if verbose: print('    Adding fgm.hmin to nc file')
        else:
            if verbose: print('fgm.hmin is None, not adding')
            
        if fgm.arrival_time is not None:        
            try:
                arrival_time = rootgrp.variables['arrival_time']
            except:
                arrival_time = rootgrp.createVariable('arrival_time','f4',('lat','lon',),
                                            fill_value=fv)
            arrival_time[:,:] = fgm.arrival_time
            arrival_time.units = 'seconds'
            if verbose: print('    Adding fgm.arrival_time to nc file')
        else:
            if verbose: print('fgm.arrival_time is None, not adding')
            
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

                
    def get_as_array(var, fgmvar=None):
        if fgmvar is None:
            fgmvar = var
        a = rootgrp.variables.get(var, None)
        if a is not None:
            if verbose: print('    Loaded %s as fgm.%s' % (var,fgmvar))
            return numpy.array(a)
        else:
            if verbose: print('    Did not find %s for fgm.%s' \
                                % (var,fgmvar))
            return None
                    
    fgm = fgmax_tools.FGmaxGrid()  # CHANGED

    with netCDF4.Dataset(fname_nc, 'r') as rootgrp:
        if verbose:
            print('Reading data to fgm from nc file',fname_nc)
            print('        nc file description: ', rootgrp.description)
            print('History:  ', rootgrp.history)


                
        x = get_as_array('lon','x')
        y = get_as_array('lat','y')
        
        if (x is None) or (y is None):
            print('*** Could not create grid')
            return None
            
        X,Y = numpy.meshgrid(x,y)
        fgm.X = X
        fgm.Y = Y
        if verbose:
            print('    Constructed fgm.X and fgm.Y')
        
        fgm.Z = get_as_array('Z')
        fgm.B = get_as_array('B')
        fgm.fgmax_point = get_as_array('fgmax_point') 
        fgm.dz = get_as_array('dz')
        fgm.h = get_as_array('h')
        fgm.s = get_as_array('s')
        fgm.hss = get_as_array('hss')
        fgm.hmin = get_as_array('hmin')
        fgm.arrival_time = get_as_array('arrival_time')
        fgm.force_dry_init = get_as_array('force_dry_init')
        
    if verbose:
        print('Returning FGmaxMaskedGrid object fgm')
    return fgm