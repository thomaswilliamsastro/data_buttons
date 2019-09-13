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
    overwrite=False,
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
        overwrite (bool, optional): Whether to create a mosaic even if 
            one already exists. Defaults to False.
        verbose (bool, optional): Print out messages during the process.
            Useful mainly for debugging purposes or large images. 
            Defaults to False.
            
    Todo:
        * Create verbose debugging statements.
    
    """

    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if filters == "all":
        filters = ["u", "g", "r", "i", "z"]

    if isinstance(filters, str):
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)

    for galaxy in galaxies:

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)

        for sdss_filter in filters:

            if os.path.exists(galaxy + "_" + sdss_filter + ".fits") and not overwrite:
                continue

            if not os.path.exists(galaxy + "/" + sdss_filter):
                os.mkdir(galaxy + "/" + sdss_filter)

            if verbose:
                print("Downloading data")

            # Montage uses its size as the length of the square, since
            # we want a radius use twice that.

            mArchiveDownload(
                "SDSS " + sdss_filter,
                galaxy,
                2 * radius.value,
                galaxy + "/" + sdss_filter,
            )

            # Mosaic all these files together.

            if verbose:
                print("Beginning mosaic")

            _ = mHdr(
                galaxy,
                2 * radius.value,
                2 * radius.value,
                galaxy + "/header.hdr",
                resolution=0.4,
            )

            tools.mosaic(
                galaxy + "/" + sdss_filter, header=galaxy + "/header.hdr", **kwargs
            )

            os.rename("mosaic/mosaic.fits", galaxy + "_" + sdss_filter + ".fits")

            # Clear out the mosaic folder.

            shutil.rmtree("mosaic/", ignore_errors=True)
            
            # Finally, convert to Jy.
            
            convert_to_jy(galaxy + "_" + sdss_filter,
                          sdss_filter)
            
def convert_to_jy(hdu_in,sdss_filter,save=True):
    
    """Convert from SDSS nanomaggies to Jy/pixel.
    
    SDSS maps are provided in convenience units of 'nanomaggies'. In 
    general, 1 nanomaggy is 3.631e-6 Jy, but there are are some offsets
    for the u and z band. This is detailed at 
    http://www.sdss3.org/dr8/algorithms/magnitudes.php.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of SDSS 
            .fits file (excluding the .fits extension), or an Astropy 
            PrimaryHDU instance (i.e. the result of ``fits.open(file)[0]``.
        sdss_filter (str): Either 'u', 'g', 'r', 'i', or 'z'.
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
    
    # Convert from nanomaggies to flux. In most bands, 1 nanomaggie = 3.631e-6 Jy.
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
    
    if save:
        fits.writeto(hdu_in+'_jy.fits',
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
