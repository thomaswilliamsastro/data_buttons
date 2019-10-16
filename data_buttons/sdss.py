# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil

import astropy.units as u
from astropy.io import fits
from MontagePy.archive import mArchiveDownload
from MontagePy.main import mHdr

from . import tools


def sdss_button(
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
    """Create an SDSS mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available SDSS data and mosaics into a final product.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        filters (str or list, optional): Any combination of 'u', 'g', 
            'r', 'i', or 'z'. If you want everything, select 'all'. 
            Defaults to 'all'.
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
        filters = ["u", "g", "r", "i", "z"]

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

        for sdss_filter in filters:
            
            if verbose:
                print('Beginning SDSS_'+sdss_filter)
                
            if not os.path.exists(galaxy + "/SDSS/"):
                os.mkdir(galaxy + "/SDSS/")

            if not os.path.exists(galaxy + "/SDSS/" + sdss_filter):
                os.mkdir(galaxy + "/SDSS/" + sdss_filter)
                
            if not os.path.exists(galaxy + "/SDSS/" + sdss_filter+'/raw'):
                os.mkdir(galaxy + "/SDSS/" + sdss_filter+'/raw')
            
            if not os.path.exists(galaxy + "/SDSS/" + sdss_filter+'/outputs'):
                os.mkdir(galaxy + "/SDSS/" + sdss_filter+'/outputs')
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
    
                # Montage uses its size as the length of the square, since
                # we want a radius use twice that.
    
                mArchiveDownload(
                    "SDSS " + sdss_filter,
                    galaxy,
                    2 * radius.value,
                    galaxy + "/SDSS/" + sdss_filter+'/raw',
                )
                
            if 2 in steps:
    
                # Mosaic all these files together.
    
                if verbose:
                    print("Beginning mosaic")
    
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/SDSS/"+sdss_filter+"/outputs/header.hdr",
                    resolution=0.4,
                    )
    
                tools.mosaic(
                    galaxy + "/SDSS/" + sdss_filter+'/raw', 
                    header=galaxy + "/SDSS/"+sdss_filter+"/outputs/header.hdr", 
                    verbose=verbose,
                    **kwargs
                    )
    
                os.rename("mosaic/mosaic.fits", 
                          galaxy + "/SDSS/" + sdss_filter + "/outputs/"+galaxy+".fits")
    
                # Clear out the mosaic folder.
    
                shutil.rmtree("mosaic/", ignore_errors=True)
                
            if 3 in steps:
                
                if verbose:
                    print('Converting to Jy')
            
                # Convert to Jy.
                
                convert_to_jy(galaxy + "/SDSS/" + sdss_filter + "/outputs/"+galaxy+".fits",
                              sdss_filter,
                              hdu_out=galaxy + "/SDSS/"+galaxy+"_" + sdss_filter+".fits")
            
def convert_to_jy(hdu_in,sdss_filter,hdu_out=None):
    
    """Convert from SDSS nanomaggies to Jy/pixel.
    
    SDSS maps are provided in convenience units of 'nanomaggies'. In 
    general, 1 nanomaggy is 3.631 x 10\ :sup:`-6`\ Jy, but there are are 
    some offsets for the u and z band. This is detailed at 
    http://www.sdss3.org/dr8/algorithms/magnitudes.php.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of SDSS 
            .fits file, or an Astropy PrimaryHDU instance (i.e. the result 
            of ``fits.open(file)[0]``).
        sdss_filter (str): Either 'u', 'g', 'r', 'i', or 'z'.
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
    
    # Convert from nanomaggies to flux. In most bands, 1 nanomaggy = 3.631e-6 Jy.
    # For the u and z band, we have magnitude corrections of +0.02 and -0.04
        
    data *= 3.631e-6
    
    if sdss_filter == 'u':
        nanomag_corr = 10**(0.04/2.51)
    elif sdss_filter == 'z':
        nanomag_corr = 10**(-0.02/2.51)
    else:
        nanomag_corr = 1
    
    data *= nanomag_corr
    
    header['BUNIT'] = 'Jy/pix'
    
    if hdu_out is not None:
        fits.writeto(hdu_out,
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
