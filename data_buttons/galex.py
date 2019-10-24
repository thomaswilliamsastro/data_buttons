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
from astroquery.ned import Ned
from MontagePy.main import mProject, mHdr

from . import tools


def galex_button(
    galaxies,
    filters="both",
    radius=None,
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
            galaxy to search for observations. Defaults to None, where
            it will query Ned to get size.
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

        obs_table = Observations.query_criteria(objectname=galaxy,
                                                radius=radius,
                                                obs_type='all',
                                                obs_collection='GALEX')
        
        # Ignore any calibration observations.
        obs_table = obs_table[obs_table['intentType'] == 'science']

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
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter+'/reprojected'):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter+'/reprojected')
                if not os.path.exists(galaxy + "/GALEX/" + galex_filter+'/weight'):
                    os.mkdir(galaxy + "/GALEX/" + galex_filter+'/weight')
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
                
            mHdr(galaxy,
                 2 * radius.value,
                 2 * radius.value,
                 galaxy + "/GALEX/"+galex_filter+"/outputs/header.hdr",
                 resolution=1.5,
                 )
                    
            if 2 in steps:

                # Read in these files and set anything more than 35
                # arcmin out to NaN.
                
                if verbose:
                    print('Performing initial weighted reprojections')

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

                    hdu.writeto(galex_file.replace('/raw/','/data/'), overwrite=True)
                    
                    # Also create a weight map (sqrt EXPTIME).
                    
                    exp_time = np.ones(hdu.data.shape)*hdu.header['EXPTIME']**0.5
                    exp_time[np.isnan(hdu.data) == True] = np.nan
                    
                    fits.writeto(galex_file.replace('/raw/','/weight/'),
                                 exp_time,hdu.header,
                                 overwrite=True)
                    
                    # And reproject each map separately using this weighting.
                    
                    mProject(galex_file.replace('/raw/','/data/'),
                             galex_file.replace('/raw/','/reprojected/'),
                             galaxy + "/GALEX/"+galex_filter+"/outputs/header.hdr",
                             weight_file=galex_file.replace('/raw/','/weight/'))

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
                    galaxy + "/GALEX/" + galex_filter+'/reprojected',
                    header=galaxy+"/GALEX/"+galex_filter+"/outputs/header.hdr",
                    verbose=verbose,
                    reproject=False,
                    haveAreas=True,
                    )
                
                os.rename(
                    "mosaic/mosaic.fits", 
                    galaxy + "/GALEX/"+galex_filter+"/outputs/"+galaxy+'.fits'
                    )

                shutil.rmtree("mosaic/", ignore_errors=True)
                    
            if 3 in steps:
                
                if verbose:
                    print('Converting to Jy')
                
                # Convert to Jy.
                
                convert_to_jy(galaxy + "/GALEX/" + galex_filter +'/outputs/'+galaxy+".fits",
                              galex_filter,
                              hdu_out=galaxy + "/GALEX/"+galaxy+'_' + galex_filter + ".fits")
                
        if original_radius is None:
            radius = None
        else:
            radius = original_radius.copy()
        
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
