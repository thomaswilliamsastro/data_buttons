# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import fnmatch
import glob
import gzip
import os
import shutil
import sys

import astropy.units as u
import numpy as np
from astropy.io import fits
from astroquery.mast import Observations
from MontagePy.main import mHdr

from . import tools


def galex_button(
    galaxies,
    filters="both",
    radius=0.2 * u.degree,
    filepath=None,
    download_data=True,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
):
    """Create a GALEX mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available GALEX data and mosaics into a final product.
    
    Because GALEX images are in counts/s and the integrations may be 
    different lengths, we convert back to a raw count, add the frames
    and convert back to counts/s at the end. This effectively weights
    the frame by exposure time.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        filters (str, optional): One of 'FUV', 'NUV', or 'both'. Selects
            which GALEX filters to create a mosaic for. Defaults to
            'both'.
        radius (astropy.units.Quantity, optional): Radius around the 
            galaxy to search for observations. Defaults to 0.2 degrees.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        download_data (bool, optional): If True, will download data from 
            MAST. Defaults to True.
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

    if filters == "both":
        filters = ["NUV", "FUV"]

    else:
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

        obs_table = Observations.query_object(galaxy, radius=radius)

        for galex_filter in filters:
            
            if verbose:
                print('Beginning GALEX '+galex_filter)
                
            if 1 in steps:
                
                # Pull out available data, and download it

                query_results = np.where(obs_table["filters"] == galex_filter)[0]
    
                # If there isn't any GALEX coverage, just skip
    
                if len(query_results) == 0:
                    print(galaxy + " missing!")
    
                    continue
                
                # We only want to print out download messages if
                # verbose is True, so redirect otherwise.
                
                if not verbose:
                    sys.stdout = open(os.devnull,'w')
                    
                dataProductsByID = Observations.get_product_list(
                    obs_table[query_results]
                    )
                  
                Observations.download_products(
                    dataProductsByID, download_dir="galex_temp/" + galaxy, mrp_only=True,
                    )
                    
                # And set back to the original for printing.
                
                if not verbose:
                    
                    sys.stdout = sys.__stdout__

                if not os.path.exists(galaxy + "/GALEX/"):
                    os.mkdir(galaxy + "/GALEX/")
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter)
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter+'/raw'):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter+'/raw')
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter+'/data'):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter+'/data')
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter+'/exptime'):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter+'/exptime')
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter+'/outputs'):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter+'/outputs')

                ext_name = {"NUV": "nd", "FUV": "fd"}[galex_filter]

                # Pull out the relevant filter files files (either *nd*
                # or *fd*), extract and move to base folder.

                matches = []
                for root, _, filenames in os.walk("galex_temp/" + galaxy):
                    for filename in fnmatch.filter(
                        filenames, "*" + ext_name + "-int.fits.gz"
                    ):
                        matches.append(os.path.join(root, filename))

                for match in matches:
                    with gzip.open(match, "rb") as f_in:
                        with open(match[:-3], "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)

                    filename = match[:-3].split('/')

                    os.rename(
                        match[:-3],
                        galaxy
                        + "/GALEX/"
                        + galex_filter
                        + "/raw/"
                        + filename[-1],
                    )
                    
                # Clean up any temporary files.

                shutil.rmtree("galex_temp/" + galaxy, ignore_errors=True)
                    
            if 2 in steps:

                # Read in these files and set anything more than 35
                # arcmin out to NaN. At the same time, convert from 
                # counts/s to raw counts so it's properly weighted in
                # the coadd later.

                galex_files = glob.glob(galaxy + "/GALEX/" + galex_filter + "/raw/*")

                for galex_file in galex_files:
                    hdu = fits.open(galex_file)[0]

                    i = np.linspace(
                        -hdu.data.shape[0] / 2, hdu.data.shape[0] / 2, hdu.data.shape[0]
                    )

                    j = np.linspace(
                        -hdu.data.shape[1] / 2, hdu.data.shape[1] / 2, hdu.data.shape[1]
                    )

                    iv, jv = np.meshgrid(i, j)
                    r = iv ** 2 + jv ** 2

                    hdu.data[r >= 1400 ** 2] = np.nan
                    
                    hdu.data *= hdu.header['EXPTIME']
                    
                    galex_filepath = galex_file.split('/raw/')
                    data_file = galex_filepath[0]+'/data/'+galex_filepath[1]

                    fits.writeto(data_file, hdu.data, hdu.header, overwrite=True)
                    
                    # Also create a map of exposure time
                    
                    exp_file = galex_filepath[0]+'/exptime/'+galex_filepath[1]
                    
                    exposure_time = np.ones(hdu.data.shape)*hdu.header['EXPTIME']
                    exposure_time[np.isnan(hdu.data) == True] = np.nan
                    fits.writeto(exp_file, exposure_time, hdu.header, overwrite=True)

                # And mosaic!

                # Montage uses its size as the length of the square,
                # since we want a radius use twice that.

                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/GALEX/"+galex_filter+"/outputs/header.hdr",
                    resolution=1.5,
                    )

                tools.mosaic(
                    galaxy + "/GALEX/" + galex_filter+'/data',
                    header=galaxy+"/GALEX/"+galex_filter+"/outputs/header.hdr",
                    verbose=verbose,
                    coadd=3
                    )
                
                os.rename(
                    "mosaic/mosaic.fits", 
                    galaxy + "/GALEX/"+galex_filter+"/outputs/"+galaxy+'_data.fits'
                    )

                tools.mosaic(
                    galaxy + "/GALEX/" + galex_filter + "/exptime",
                    header=galaxy+"/GALEX/"+galex_filter+"/outputs/header.hdr",
                    verbose=verbose,
                    coadd=3
                    )

                os.rename(
                    "mosaic/mosaic.fits", 
                    galaxy + "/GALEX/" + galex_filter +'/outputs/'+galaxy+"_exptime.fits"
                    )

                shutil.rmtree("mosaic/", ignore_errors=True)
                    
                # Convert back to counts/sec
                
                hdu = fits.open(galaxy + "/GALEX/" + galex_filter +'/outputs/'+galaxy+"_data.fits")[0]
                hdu_exp = fits.open(galaxy + "/GALEX/" + galex_filter +'/outputs/'+galaxy+"_exptime.fits")[0]
                
                hdu.data /= hdu_exp.data
                
                fits.writeto(galaxy + "/GALEX/" + galex_filter +'/outputs/'+galaxy+".fits",
                             hdu.data,hdu.header,
                             overwrite=True)
                    
            if 3 in steps:
                
                if verbose:
                    print('Converting to Jy')
                
                # Convert to Jy.
                
                convert_to_jy(galaxy + "/GALEX/" + galex_filter +'/outputs/'+galaxy+".fits",
                              galex_filter,
                              hdu_out=galaxy + "/GALEX/"+galaxy+'_' + galex_filter + ".fits")
        
def convert_to_jy(hdu_in,galex_filter,hdu_out=None):
    
    """Convert from GALEX counts/s to Jy/pixel.
    
    GALEX maps are provided in units of counts/s. These can be converted
    to more helpful units via the conversion factors given at
    https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of GALEX 
            .fits file, or an Astropy PrimaryHDU instance (i.e. the result 
            of ``fits.open(file)[0]``).
        galex_filter (str): Either 'NUV' or 'FUV'.
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
    
    # Conversion factors here from Clark+ (2017)
    
    conversion_factor = {'FUV':1.076e-4,
                         'NUV':3.373e-5}[galex_filter]
    
    data *= conversion_factor
    header['BUNIT'] = 'Jy/pix'
    
    if hdu_out is not None:
        fits.writeto(hdu_out,
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
