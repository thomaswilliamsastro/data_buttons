# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil

import astropy.units as u
from MontagePy.archive import mArchiveDownload
from MontagePy.main import mHdr

from . import tools


def sdss_button(
        galaxies,
        filters='all',
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
    
    """

    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if filters == 'all':
        filters = ['u', 'g', 'r', 'i', 'z']

    if isinstance(filters, str):
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)

    for galaxy in galaxies:

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)

        for sdss_filter in filters:

            if os.path.exists(galaxy + '_' + sdss_filter + '.fits') and \
                    not overwrite:
                continue

            if not os.path.exists(galaxy + '/' + sdss_filter):
                os.mkdir(galaxy + '/' + sdss_filter)

            if verbose:
                print('Downloading data')
                
            # Montage uses its size as the length of the square, since 
            # we want a radius use twice that.

            mArchiveDownload(
                'SDSS ' + sdss_filter,
                galaxy, 2*radius.value,
                galaxy + '/' + sdss_filter,
                )

            # Mosaic all these files together.

            if verbose:
                print('Beginning mosaic')
                
            _ = mHdr(galaxy,
                     2*radius.value, 
                     2*radius.value,
                     galaxy+'/header.hdr',
                     resolution=0.4,
                     )

            tools.mosaic(galaxy + '/' + sdss_filter,
                         header=galaxy+'/header.hdr',
                         **kwargs
                         )

            os.rename('mosaic/mosaic.fits',
                      galaxy + '_' + sdss_filter + '.fits')

            # Clear out the mosaic folder.

            shutil.rmtree('mosaic/',
                          ignore_errors=True,
                          )
