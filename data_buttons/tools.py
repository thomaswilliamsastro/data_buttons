# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil
import warnings

import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from MontagePy.archive import *
from MontagePy.main import *
from photutils import make_source_mask


def mosaic(input_folder, header=None, output_folder="mosaic", background_match=True,
           verbose=False,**kwargs):

    """Mosaic together a folder full of .fits files.
    
    Takes a folder and runs through the number of steps that MontagePy
    requires to combine all the .fits files within that folder. Also,
    optionally performs background modelling and matching (which most
    users will want to set to true). Note that this is not a background
    subtraction!
    
    Args:
        input_folder (str): Folder of raw files to mosaic.
        header (str, optional): Output from mHdr. If not specified, 
            will mosaic all of the images together no matter the
            overlap.
        output_folder (str, optional): Working folder for mosaicking.
            Defaults to 'mosaic'.
        background_match (bool, optional): Whether to perform background
            matching steps while mosaicking. Defaults to True.
        verbose (bool, optional): Print out verbose statements during
            the mosaicking process. Useful for debugging. Defaults to
            False.
            
    """

    # Create the working folders we'll need

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    os.mkdir(output_folder + "/projected")
    if background_match:
        os.mkdir(output_folder + "/diffs")
        os.mkdir(output_folder + "/corrected")

    # Make an optimum header for these images

    mImgtbl(input_folder, output_folder + "/images.tbl")

    if header is None:
        mMakeHdr(output_folder + "/images.tbl", output_folder + "/header.hdr")
    else:
        shutil.copy(header, output_folder + "/header.hdr")

    # Project the original images to this header
    
    if verbose:
        print('Reprojecting images to optimum header')

    mProjExec(
        input_folder,
        output_folder + "/images.tbl",
        output_folder + "/header.hdr",
        projdir=output_folder + "/projected",
        quickMode=False,
    )

    mImgtbl(output_folder + "/projected", output_folder + "/images.tbl")

    # If selected, perform background matching.

    if background_match:
        
        if verbose:
            print('Calculating overlaps')
        
        mOverlaps(output_folder + "/images.tbl", output_folder + "/diffs.tbl")
        
        if verbose:
            print('Fitting overlap differences')
        
        mDiffFitExec(
            output_folder + "/projected",
            output_folder + "/diffs.tbl",
            output_folder + "/header.hdr",
            output_folder + "/diffs",
            output_folder + "/fits.tbl",
        )
        
        if verbose:
            print('Calulating corrections')
        
        mBgModel(
            output_folder + "/images.tbl",
            output_folder + "/fits.tbl",
            output_folder + "/corrections.tbl",
        )
        
        if verbose:
            print('Matching backgrounds')
        
        mBgExec(
            output_folder + "/projected",
            output_folder + "/images.tbl",
            output_folder + "/corrections.tbl",
            output_folder + "/corrected",
        )
        mImgtbl(output_folder + "/corrected", output_folder + "/images.tbl")

    # Finally, coadd the images

    if not background_match:
        folder = output_folder + "/projected"
    else:
        folder = output_folder + "/corrected"
        
    if verbose:
        print('Co-adding images')

    mAdd(
        folder,
        output_folder + "/images.tbl",
        output_folder + "/header.hdr",
        output_folder + "/mosaic.fits",
        **kwargs
    )

    # Remove the temp folders we've made along the way

#     shutil.rmtree(output_folder + "/projected", ignore_errors=True)
#     if background_match:
#         shutil.rmtree(output_folder + "/diffs", ignore_errors=True)
#         shutil.rmtree(output_folder + "/corrected", ignore_errors=True)
        
def calculate_background_median(data,sigma=3,npixels=5,maxiters=20,
                                **kwargs):
    
    """Calculate background median for data.
    
    Calculates the background median for an image by iteratively sigma
    clipping. It first creates a source mask using segmentation and binary
    dilation (see https://photutils.readthedocs.io/en/stable/api/photutils.make_source_mask.html)
    for more details, before iteratively calculating the sigma-clipped 
    median of the data.
    
    Args:
        data (str or numpy.ndarray or astropy.io.fits.PrimaryHDU): Data 
            to calculate the background median for. If provided as a 
            string, this is interpreted as a filename for a .fits HDU.
        sigma (float, optional): Level to perform sigma-clipping to. 
            Defaults to 3.
        npixels (int, optional): Number of connected pixels greater than 
            sigma to consider for a pixel to be part of a source when 
            masking. Defaults to 5.
        maxiters (int,optional): The maximum number of sigma-clipping 
            iterations to perform. Defaults to 20.
            
    Returns:
        float: Background median of the data.
    
    """
    
    if isinstance(data,str):
        hdu = fits.open(data)[0]
        data = hdu.data 
    elif isinstance(data,fits.PrimaryHDU):
        data = data.data
    
    mask = make_source_mask(data,snr=sigma,npixels=npixels,
                            **kwargs)
    _,median,_ = sigma_clipped_stats(data,mask=mask,sigma=sigma,
                                     maxiters=maxiters,**kwargs)
    
    return median 

def model_background(data,poly_order=5,sigma=2,npixels=5,
                     **kwargs):
    
    """Model the background as a 2D polynomial.
    
    After masking sources, fits a 2D polynomial to the remaining data.
    This is useful for images where the background ripples over small
    angular scales, such as 2MASS or SDSS data.
    
    Args:
        data (str or numpy.ndarray or astropy.io.fits.PrimaryHDU): Data 
            to calculate the background median for. If provided as a 
            string, this is interpreted as a filename for a .fits HDU.
        poly_order (int, optional): The 2D polynomial order to fit the
            background to. Defaults to 5.
        sigma (float, optional): Level to set mask minimum SNR to. 
            Defaults to 2.
        npixels (int, optional): Number of connected pixels greater than 
            sigma to consider for a pixel to be part of a source when 
            masking. Defaults to 5.
    
    Returns:
        numpy.ndarray: The modelled background; same size as the input
            data.
    
    """
    
    if isinstance(data,str):
        hdu = fits.open(data)[0]
        data = hdu.data 
    elif isinstance(data,fits.PrimaryHDU):
        data = data.data
    
    y,x = np.mgrid[:data.shape[0],:data.shape[1]]
    
    # Create a mask for the data.
    
    mask = make_source_mask(data,snr=sigma,npixels=npixels,
                            **kwargs)
    data_masked = data.copy()
    data_masked[mask == True] = np.nan
    
    # Remove any NaNs (i.e. data missing or data masked).
    
    not_nan = np.where( np.isnan(data_masked) == False )
    
    x_fit = x[not_nan]
    y_fit = y[not_nan]
    data_fit = data[not_nan]
    
    # Fit background to the data.

    p_init = models.Polynomial2D(poly_order)
    fit_p = fitting.LevMarLSQFitter()
    
    with warnings.catch_warnings():
        
        # Ignore model linearity warning from the fitter.
        warnings.simplefilter('ignore')
        p = fit_p(p_init, x_fit, y_fit, data_fit)
    
    background = p(x,y)
    
    return background
