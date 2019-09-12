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
            not os.path.exists(galaxy + "_FUV.fits")
            or not os.path.exists(galaxy + "_NUV.fits")
            or overwrite
        ):

            obs_table = Observations.query_object(galaxy, radius=radius)

            # Pull out available data, and download it

            query_results = np.where(obs_table["filters"] == filters[0])[0]

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

            for galex_filter in filters:

                if (
                    os.path.exists(galaxy + "_" + galex_filter + ".fits")
                    and not overwrite
                ):
                    continue

                if not os.path.exists(galaxy + "/" + galex_filter):
                    os.mkdir(galaxy + "/" + galex_filter)

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
                        + "/"
                        + galaxy
                        + "_"
                        + str(filename_ext)
                        + ".fits",
                    )

                    filename_ext += 1

                # Read in these files and set anything more than 35
                # arcmin out to NaN

                galex_files = glob.glob(galaxy + "/" + galex_filter + "/*")

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

                    fits.writeto(galex_file, hdu.data, hdu.header, overwrite=True)

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
                        galaxy + "/" + galex_filter,
                        header=galaxy + "/header.hdr",
                        **kwargs
                    )

                    # Move the mosaic out and rename

                    os.rename(
                        "mosaic/mosaic.fits", galaxy + "_" + galex_filter + ".fits"
                    )

                    shutil.rmtree("mosaic/", ignore_errors=True)

                else:

                    os.rename(
                        galaxy + "/" + galex_filter + "/" + galaxy + "_0.fits",
                        galaxy + "_" + galex_filter + ".fits",
                    )

        # Clean up any temporary files.

        shutil.rmtree("galex_temp/" + galaxy, ignore_errors=True)
