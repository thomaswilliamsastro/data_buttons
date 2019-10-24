# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil
import numpy as np

import astropy.units as u
from astropy.io import fits
from astroquery.ned import Ned
from MontagePy.archive import mArchiveDownload
from MontagePy.main import mHdr

from . import tools

# New imports


def wise_button(
    galaxies,
    filters="all",
    radius=None,
    filepath=None,
    download_data=True,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
    **kwargs
):
    
    """Create a WISE mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available WISE data and mosaics into a final product.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        filters (str or list, optional): Any combination of '1', '2', 
            '3', and '4'. If you want everything, select 'all'. Defaults 
            to 'all'.
        radius (astropy.units.Quantity, optional): Radius around the 
            galaxy to search for observations. Defaults to None, where
            it will query Ned to get size.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        download_data (bool, optional): If True, will download data using
            MontagePy. Defaults to True.
        create_mosaic (bool, optional): Switching this to True will 
            mosaic data as appropriate. Defaults to True.
        jy_conversion (bool, optional): Convert the mosaicked file from
            raw units to Jy/pix. Defaults to True.
        verbose (bool, optional): Print out messages during the process.
            Useful mainly for debugging purposes or large images. 
            Defaults to False.
    
    """
    
    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if filters == "all":
        filters = ['1','2','3','4']

    if isinstance(filters, str):
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)
        
    if radius is not None:
        original_radius = radius.copy()
    else:
        original_radius = None
        
    steps = []
    
    if download_data:
        steps.append(1)
    if create_mosaic:
        steps.append(2)
    if jy_conversion:
        steps.append(3)

    for galaxy in galaxies:
        
        if verbose:
            print('Beginning '+galaxy)
            
        if radius is None:
        
            try:
 
                size_query = Ned.get_table(galaxy,table='diameters')
                radius = 1.2*np.max(size_query['NED Major Axis'])/2*u.arcsec
                radius = radius.to(u.deg)
     
            except:
                
                raise Warning(galaxy+' not resolved by Ned, using 0.2deg radius.')
                radius = 0.2*u.degree

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)

        for wise_filter in filters:
            
            if verbose:
                print('Beginning W'+wise_filter)
            
            if not os.path.exists(galaxy + "/WISE"):
                os.mkdir(galaxy + "/WISE")
                
            if not os.path.exists(galaxy + "/WISE/W" + wise_filter):
                os.mkdir(galaxy + "/WISE/W" + wise_filter)
                
            if not os.path.exists(galaxy + "/WISE/W" + wise_filter+'/raw'):
                os.mkdir(galaxy + "/WISE/W" + wise_filter+'/raw')
                
            if not os.path.exists(galaxy + "/WISE/W" + wise_filter+'/outputs'):
                os.mkdir(galaxy + "/WISE/W" + wise_filter+'/outputs')
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
    
                # Montage uses its size as the length of the square, 
                # since we want a radius use twice that.
    
                mArchiveDownload(
                    "WISE " + wise_filter,
                    galaxy,
                    2 * radius.value,
                    galaxy + "/WISE/W" + wise_filter+'/raw',
                )
                
            if 2 in steps:
                
                # Mosaic all these files together.
    
                if verbose:
                    print("Beginning mosaic")
    
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/WISE/W" + wise_filter+"/outputs/header.hdr",
                    resolution=1.375,
                    )
    
                tools.mosaic(
                    galaxy + "/WISE/W" + wise_filter+'/raw', 
                    header=galaxy + "/WISE/W" + wise_filter+"/outputs/header.hdr", 
                    verbose=verbose,
                    **kwargs
                    )
    
                os.rename("mosaic/mosaic.fits", 
                          galaxy + "/WISE/W" + wise_filter + "/outputs/"+galaxy+".fits")
    
                # Clear out the mosaic folder.
    
                shutil.rmtree("mosaic/", ignore_errors=True)
            
            if 3 in steps:
                
                if verbose:
                    print('Converting to Jy')
                    
                # Convert to Jy.
            
                convert_to_jy(galaxy + "/WISE/W" + wise_filter + "/outputs/"+galaxy+".fits",
                              wise_filter,
                              hdu_out=galaxy + "/WISE/"+galaxy+"_W" + wise_filter+".fits")
            
        if original_radius is None:
            radius = None
        else:
            radius = original_radius.copy()
            
def convert_to_jy(hdu_in,wise_filter,hdu_out=None):
    
    """Convert from WISE DN to Jy/pixel.
    
    WISE maps are provided in convenience units of data numbers (DN). The
    constants to convert them are given in Table 1 of
    http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec2_3f.html.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of WISE 
            .fits file, or an Astropy PrimaryHDU instance (i.e. the
            result of ``fits.open(file)[0]``).
        wise_filter (str): Either '1', '2', '3', or '4'.
        hdu_out (str, optional): If not None, will save the converted HDU
            out with this filename. Defaults to None.
        
    Returns:
        astropy.io.fits.PrimaryHDU: The HDU in units of Jy/pix.
    
    """
    
    if isinstance(hdu_in,str):
        hdu = fits.open(hdu_in)[0]
    else:
        hdu = hdu_in.copy()
    
    data = hdu.data.copy()
    header = hdu.header.copy()
    
    # Convert from DN to Jy
        
    dn_factor = {'1':1.935e-6,
                 '2':2.7048e-6,
                 '3':2.9045e-6,
                 '4':5.2269e-5}[wise_filter]
    
    data *= dn_factor
    
    header['BUNIT'] = 'Jy/pix'
    
    if hdu_out is not None:
        fits.writeto(hdu_out,
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
            
            