# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil

import astropy.units as u
from astropy.io import fits
from MontagePy.archive import mArchiveDownload
from MontagePy.main import mHdr

from . import tools

# New imports


def wise_button(
    galaxies,
    filters="all",
    radius=0.2 * u.degree,
    filepath=None,
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
            galaxy to search for observations. Defaults to 0.2 degrees.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        create_mosaic (bool, optional): Switching this to True will 
            download data and mosaic as appropriate. You may wish to set
            this to False if you've already downloaded the data previously.
            Defaults to True.
        jy_conversion (bool, optional): Convert the mosaicked file from
            raw units to Jy/pix. Defaults to True.
        verbose (bool, optional): Print out messages during the process.
            Useful mainly for debugging purposes or large images. 
            Defaults to False.
            
    Todo:
    
    """
    
    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if filters == "all":
        filters = ['1','2','3','4']

    if isinstance(filters, str):
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)
        
    steps = []
    
    if create_mosaic:
        steps.append(1)
    if jy_conversion:
        steps.append(2)

    for galaxy in galaxies:
        
        if verbose:
            print('Beginning '+galaxy)

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)

        for wise_filter in filters:
            
            if verbose:
                print('Beginning W'+wise_filter)
            
            if not os.path.exists(galaxy + "/W" + wise_filter):
                os.mkdir(galaxy + "/W" + wise_filter)
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
    
                # Montage uses its size as the length of the square, 
                # since we want a radius use twice that.
    
                mArchiveDownload(
                    "WISE " + wise_filter,
                    galaxy,
                    2 * radius.value,
                    galaxy + "/W" + wise_filter,
                )
                
                # Mosaic all these files together.
    
                if verbose:
                    print("Beginning mosaic")
    
                _ = mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/header.hdr",
                    resolution=1.375,
                )
    
                tools.mosaic(
                    galaxy + "/W" + wise_filter, 
                    header=galaxy + "/header.hdr", 
                    **kwargs
                )
    
                os.rename("mosaic/mosaic.fits", 
                          galaxy + "_W" + wise_filter + ".fits")
    
                # Clear out the mosaic folder.
    
                shutil.rmtree("mosaic/", ignore_errors=True)
            
            # Convert to Jy.
            
            if 2 in steps:
                
                print('Converting to Jy')
            
                convert_to_jy(galaxy + "_W" + wise_filter,
                              wise_filter)
            
def convert_to_jy(hdu_in,wise_filter,save=True):
    
    """Convert from WISE DN to Jy/pixel.
    
    WISE maps are provided in convenience units of data numbers (DN). The
    constants to convert them are given in Table 1 of
    http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec2_3f.html.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of WISE 
            .fits file (excluding the .fits extension), or an Astropy 
            PrimaryHDU instance (i.e. the result of ``fits.open(file)[0]``).
        galex_filter (str): Either '1', '2', '3', or '4'.
        save (bool, optional): Save out the converted file. It'll save
            the original file with an appended '_jy'. Defaults to True.
        
    Returns:
        hdu_pixel: The HDU in units of Jy/pix.
    
    """
    
    if isinstance(hdu_in,str):
        hdu = fits.open(hdu_in+'.fits')[0]
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
    
    if save:
        fits.writeto(hdu_in+'_jy.fits',
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
            
            