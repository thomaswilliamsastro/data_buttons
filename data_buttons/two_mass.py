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
    download_data=True,
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
        filters = ['J', 'H', 'K']

    if isinstance(filters, str):
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)
        
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

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)

        for two_mass_filter in filters:
            
            if verbose:
                print('Beginning 2MASS '+two_mass_filter)
                
            if not os.path.exists(galaxy + "/2MASS"):
                os.mkdir(galaxy + "/2MASS")
            
            if not os.path.exists(galaxy + "/2MASS/" + two_mass_filter):
                os.mkdir(galaxy + "/2MASS/" + two_mass_filter)
                
            if not os.path.exists(galaxy + "/2MASS/" + two_mass_filter+"/raw"):
                os.mkdir(galaxy + "/2MASS/" + two_mass_filter+"/raw")
                
            if not os.path.exists(galaxy + "/2MASS/" + two_mass_filter+"/data"):
                os.mkdir(galaxy + "/2MASS/" + two_mass_filter+"/data")
                
            if not os.path.exists(galaxy + "/2MASS/" + two_mass_filter+"/outputs"):
                os.mkdir(galaxy + "/2MASS/" + two_mass_filter+"/outputs")
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
    
                # Montage uses its size as the length of the square, 
                # since we want a radius use twice that.
    
                mArchiveDownload(
                    "2MASS " + two_mass_filter,
                    galaxy,
                    2 * radius.value,
                    galaxy + "/2MASS/" + two_mass_filter+"/raw",
                )
                
                # We need to now convert these into magnitudes.
                
                two_mass_files = glob.glob(galaxy+'/2MASS/'+ two_mass_filter+'/raw/*')
                
                for two_mass_file in two_mass_files:
                    
                    hdu = fits.open(two_mass_file)[0]
                    magzp = hdu.header['MAGZP']
                    
                    hdu.data = magzp - 2.5*np.log10(hdu.data)
                    
                    fits.writeto(two_mass_file.replace('/raw/','/data/'),
                                 hdu.data,hdu.header,
                                 overwrite=True)
                    
            if 2 in steps:
                
                # Mosaic all these files together.
    
                if verbose:
                    print("Beginning mosaic")
    
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/2MASS/"+ two_mass_filter+"/outputs/header.hdr",
                    resolution=1,
                    )
    
                tools.mosaic(
                    galaxy + "/2MASS/" + two_mass_filter+"/data", 
                    header=galaxy + "/2MASS/"+ two_mass_filter+"/outputs/header.hdr", 
                    verbose=verbose,
                    **kwargs
                    )
    
                os.rename("mosaic/mosaic.fits", 
                          galaxy + "/2MASS/"+ two_mass_filter+"/outputs/"+galaxy+ ".fits")
    
                # Clear out the mosaic folder.
    
                shutil.rmtree("mosaic/", ignore_errors=True)
            
            if 3 in steps:
                
                if verbose:
                    print('Converting to Jy')
                    
                # Convert to Jy.
            
                convert_to_jy(galaxy + "/2MASS/"+ two_mass_filter+"/outputs/"+galaxy+ ".fits",
                              two_mass_filter,
                              galaxy + "/2MASS/"+galaxy+"_"+ two_mass_filter+".fits")
            
def convert_to_jy(hdu_in,two_mass_filter,hdu_out=None):
    
    """Convert from Vega magnitudes to Jy/pixel.
    
    2MASS maps are provided in convenience units of data numbers (DN).
    During the mosaicking process, these are converted to Vega mags
    since the magnitude zero point varies between frames. The constants 
    for converting magnitudes to Jy are given at
    https://old.ipac.caltech.edu/2mass/releases/allsky/faq.html#jansky.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of 2MASS 
            .fits file, or an Astropy PrimaryHDU instance (i.e. the result 
            of ``fits.open(file)[0]``).
        two_mass_filter (str): Either 'J', 'H', or 'K'.
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
    
    # Convert from mag to Jy
        
    f_0 = {'J':1594,
           'H':1024,
           'K':666.7}[two_mass_filter]
    
    data = f_0*10**(-0.4*data)
    
    header['BUNIT'] = 'Jy/pix'
    
    if hdu_out is not None:
        fits.writeto(hdu_out,
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
            
            