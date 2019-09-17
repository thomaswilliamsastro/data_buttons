# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil
import glob

import astropy.units as u
from astropy.io import fits
import numpy as np
from MontagePy.archive import mArchiveDownload
from MontagePy.main import mHdr

from . import tools

# New imports


def two_mass_button(
    galaxies,
    filters="all",
    radius=0.2 * u.degree,
    filepath=None,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
    **kwargs
):
    
    """Create a 2MASS mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available 2MASS data and mosaics into a final product.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        filters (str or list, optional): Any combination of 'J', 'H', 
            and 'K'. If you want everything, select 'all'. Defaults 
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
        filters = ['J', 'H', 'K']

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

        for two_mass_filter in filters:
            
            if verbose:
                print('Beginning 2MASS '+two_mass_filter)
            
            if not os.path.exists(galaxy + "/2MASS_" + two_mass_filter):
                os.mkdir(galaxy + "/2MASS_" + two_mass_filter)
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
    
                # Montage uses its size as the length of the square, 
                # since we want a radius use twice that.
    
                mArchiveDownload(
                    "2MASS " + two_mass_filter,
                    galaxy,
                    2 * radius.value,
                    galaxy + "/2MASS_" + two_mass_filter,
                )
                
                # We need to now convert these into magnitudes.
                
                two_mass_files = glob.glob(galaxy+'/2MASS_'+ two_mass_filter+'/*')
                
                for two_mass_file in two_mass_files:
                    
                    hdu = fits.open(two_mass_file)[0]
                    magzp = hdu.header['MAGZP']
                    
                    hdu.data = magzp - 2.5*np.log10(hdu.data)
                    
                    fits.writeto(two_mass_file,
                                 hdu.data,hdu.header,
                                 overwrite=True)
                
                # Mosaic all these files together.
    
                if verbose:
                    print("Beginning mosaic")
    
                _ = mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/header.hdr",
                    resolution=1,
                )
    
                tools.mosaic(
                    galaxy + "/2MASS_" + two_mass_filter, 
                    header=galaxy + "/header.hdr", 
                    **kwargs
                )
    
                os.rename("mosaic/mosaic.fits", 
                          galaxy + "_2MASS_" + two_mass_filter + ".fits")
    
                # Clear out the mosaic folder.
    
                shutil.rmtree("mosaic/", ignore_errors=True)
            
            if 2 in steps:
                
                if verbose:
                    print('Converting to Jy')
                    
                # Convert to Jy.
            
                convert_to_jy(galaxy + "_2MASS_" + two_mass_filter,
                              two_mass_filter)
            
def convert_to_jy(hdu_in,two_mass_filter,save=True):
    
    """Convert from Vega magnitudes to Jy/pixel.
    
    2MASS maps are provided in convenience units of data numbers (DN).
    During the mosaicking process, these are converted to Vega mags
    since the magnitude zero point varies between frames. The constants 
    for converting magnitudes to Jy are given at
    https://old.ipac.caltech.edu/2mass/releases/allsky/faq.html#jansky.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of 2MASS 
            .fits file (excluding the .fits extension), or an Astropy 
            PrimaryHDU instance (i.e. the result of ``fits.open(file)[0]``).
        two_mass_filter (str): Either 'J', 'H', or 'K'.
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
    
    # Convert from mag to Jy
        
    f_0 = {'J':1594,
           'H':1024,
           'K':666.7}[two_mass_filter]
    
    data = f_0*10**(-0.4*data)
    
    header['BUNIT'] = 'Jy/pix'
    
    if save:
        fits.writeto(hdu_in+'_jy.fits',
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
            
            