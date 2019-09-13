# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import fnmatch
import glob
import gzip
import os
import shutil

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
    overwrite=False,
    **kwargs
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
            galaxy to search for observations. Defaults to 0.2 degrees,
            the same as Astroquery.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        overwrite (bool, optional): Whether to create a mosaic even if 
            one already exists. Defaults to False.
            
    Todo:
        * Create verbose debugging statements.
        * Find a way to turn off downloading files information.
    
    """

    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if filters == "both":
        filters = ["NUV", "FUV"]

    else:
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)

    for galaxy in galaxies:

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)

        if (
            not os.path.exists(galaxy + "_FUV_jy.fits")
            or not os.path.exists(galaxy + "_NUV_jy.fits")
            or overwrite
        ):

            obs_table = Observations.query_object(galaxy, radius=radius)

            for galex_filter in filters:
                
                # Pull out available data, and download it

                query_results = np.where(obs_table["filters"] == galex_filter)[0]
    
                # If there isn't any GALEX coverage, just skip
    
                if len(query_results) == 0:
                    print(galaxy + " missing!")
    
                    continue
    
                for query_result in query_results:
                    dataProductsByID = Observations.get_product_list(
                        obs_table[query_result]
                    )
                    _ = Observations.download_products(
                        dataProductsByID, download_dir="galex_temp/" + galaxy, mrp_only=True
                    )

                if (
                    os.path.exists(galaxy + "_" + galex_filter + ".fits")
                    and not overwrite
                ):
                    continue

                if not os.path.exists(galaxy + "/" + galex_filter):
                    os.mkdir(galaxy + "/" + galex_filter)
                if not os.path.exists(galaxy + "/" + galex_filter+'/int'):
                    os.mkdir(galaxy + "/" + galex_filter+'/int')
                if not os.path.exists(galaxy + "/" + galex_filter+'/exp'):
                    os.mkdir(galaxy + "/" + galex_filter+'/exp')

                ext_name = {"NUV": "nd", "FUV": "fd"}[galex_filter]

                # Pull out the relevant filter files files (either *nd*
                # or *fd*), extract and move to base folder.

                filename_ext = 0

                matches = []
                for root, _, filenames in os.walk("galex_temp/" + galaxy):
                    for filename in fnmatch.filter(
                        filenames, "*" + ext_name + "-int.fits.gz"
                    ):
                        matches.append(os.path.join(root, filename))

                filename_ext = 0

                for match in matches:
                    with gzip.open(match, "rb") as f_in:
                        with open(match[:-3], "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)

                    os.rename(
                        match[:-3],
                        galaxy
                        + "/"
                        + galex_filter
                        + "/int/"
                        + galaxy
                        + "_"
                        + str(filename_ext)
                        + ".fits",
                    )

                    filename_ext += 1

                # Read in these files and set anything more than 35
                # arcmin out to NaN. At the same time, convert from 
                # counts/s to raw counts so it's properly weighted in
                # the coadd later.

                galex_files = glob.glob(galaxy + "/" + galex_filter + "/int/*")

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

                    fits.writeto(galex_file, hdu.data, hdu.header, overwrite=True)
                    
                    # Also create a map of exposure time
                    
                    exp_file = galex_file.split('/int/')
                    exp_file = exp_file[0]+'/exp/'+exp_file[1]
                    
                    exposure_time = np.ones(hdu.data.shape)*hdu.header['EXPTIME']
                    exposure_time[np.isnan(hdu.data) == True] = np.nan
                    fits.writeto(exp_file, exposure_time, hdu.header, overwrite=True)

                # And mosaic!

                if len(galex_files) > 1:

                    # Montage uses its size as the length of the square,
                    # since we want a radius use twice that.

                    _ = mHdr(
                        galaxy,
                        2 * radius.value,
                        2 * radius.value,
                        galaxy + "/header.hdr",
                        resolution=1.5,
                    )

                    tools.mosaic(
                        galaxy + "/" + galex_filter+'/int',
                        header=galaxy + "/header.hdr",
                        coadd=3
                    )
                    
                    os.rename(
                        "mosaic/mosaic.fits", galaxy + "_" + galex_filter + ".fits"
                    )

                    tools.mosaic(
                        galaxy + "/" + galex_filter + "/exp",
                        header=galaxy + "/header.hdr",
                        coadd=3
                    )

                    os.rename(
                        "mosaic/mosaic.fits", galaxy + "_" + galex_filter + "_exp.fits"
                    )

                    shutil.rmtree("mosaic/", ignore_errors=True)

                else:

                    os.rename(
                        galaxy + "/" + galex_filter + "/" + galaxy + "/exp/_0.fits",
                        galaxy + "_" + galex_filter + ".fits",
                    )
                    os.rename(
                        galaxy + "/" + galex_filter + "/" + galaxy + "/int/_0.fits",
                        galaxy + "_" + galex_filter + "_exp.fits",
                    )
                    
                # Convert back to counts/sec
                
                hdu = fits.open(galaxy + "_" + galex_filter + ".fits")[0]
                hdu_exp = fits.open(galaxy + "_" + galex_filter + "_exp.fits")[0]
                
                hdu.data /= hdu_exp.data
                
                fits.writeto(galaxy + "_" + galex_filter + ".fits",
                             hdu.data,hdu.header,
                             overwrite=True)
                    
                # Convert to Jy.
                
                convert_to_jy(galaxy + "_" + galex_filter,
                              galex_filter)

                # Clean up any temporary files.
        
                shutil.rmtree("galex_temp/" + galaxy, ignore_errors=True)
        
def convert_to_jy(hdu_in,galex_filter,save=True):
    
    """Convert from GALEX counts/s to Jy/pixel.
    
    GALEX maps are provided in units of counts/s. These can be converted
    to more helpful units via the conversion factors given at
    https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of GALEX 
            .fits file (excluding the .fits extension), or an Astropy 
            PrimaryHDU instance (i.e. the result of ``fits.open(file)[0]``.
        galex_filter (str): Either 'NUV' or 'FUV'.
        save (bool, optional): Save out the converted file. It'll save
            the original file with an appended '_jy'. Defaults to True.
        
    Returns:
        astropy.io.fits.PrimaryHDU: The HDU in units of Jy/pix.
    
    """
    
    if isinstance(hdu_in,str):
        hdu = fits.open(hdu_in+'.fits')[0]
    else:
        hdu = hdu_in.copy()
    
    data = hdu.data.copy()
    header = hdu.header.copy()
    
    # Conversion factors here from Clark+ (2017)
    
    conversion_factor = {'FUV':1.076e-4,
                         'NUV':3.373e-5}[galex_filter]
    
    data *= conversion_factor
    header['BUNIT'] = 'Jy/pix'
    
    if save:
        fits.writeto(hdu_in+'_jy.fits',
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
