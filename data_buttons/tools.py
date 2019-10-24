# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil
import warnings
import glob
import wget

import numpy as np
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.nddata.utils import Cutout2D
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astropy.wcs import WCS
from MontagePy.archive import *
from MontagePy.main import *
from photutils import make_source_mask


def mosaic(input_folder, header=None, output_folder="mosaic", reproject=True,
           background_match=True, haveAreas=False,verbose=False,):

    """Mosaic together a folder full of .fits files.
    
    Takes a folder and runs through the number of steps that MontagePy
    requires to combine all the .fits files within that folder. Also,
    optionally performs background modelling and matching (which most
    users will want to set to true). Note that this is not a background
    subtraction! Can also be set to weight the final added images, for
    instance if you have input files of different exposure times.
    
    Args:
        input_folder (str): Folder of raw files to mosaic.
        header (str, optional): Output from mHdr. If not specified, 
            will mosaic all of the images together no matter the
            overlap.
        output_folder (str, optional): Working folder for mosaicking.
            Defaults to 'mosaic'.
        reproject (bool, optional): Whether to reproject input images or
            not. If already on same pixel scale, you can skip this. In 
            this case, it will move files from input_folder to the projected
            folder. Defaults to True.
        background_match (bool, optional): Whether to perform background
            matching steps while mosaicking. Defaults to True.
        haveAreas (bool, optional): If weighting a mosaic, you can use
            _area.fits files to effectively weight properly. If true, 
            before co-adding the area files will be moved alongside
            the data files. Defaults to False.
        verbose (bool, optional): Print out verbose statements during
            the mosaicking process. Useful for debugging. Defaults to
            False.
            
    """

    # Create the working folders we'll need

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
 
    if not os.path.exists(output_folder + "/projected") and reproject:
        os.mkdir(output_folder + "/projected")
         
    if background_match:
        if not os.path.exists(output_folder + "/diffs"):
            os.mkdir(output_folder + "/diffs")
        if not os.path.exists(output_folder + "/corrected"):
            os.mkdir(output_folder + "/corrected")

    # Make an optimum header for these images
    
    mImgtbl(input_folder, output_folder + "/images.tbl",)

    if header is None:
        mMakeHdr(output_folder + "/images.tbl", output_folder + "/header.hdr")
    else:
        shutil.copy(header, output_folder + "/header.hdr")

    # Project the original images to this header
    
    if reproject:
        
        if verbose:
            print('Reprojecting images to optimum header')
    
        mProjExec(
            input_folder,
            output_folder + "/images.tbl",
            output_folder + "/header.hdr",
            projdir=output_folder + "/projected",
            quickMode=False,
        )

    else:
        
        if verbose:
            print('Moving existing reprojected images')
            
        shutil.copytree(input_folder,
                        output_folder+"/projected")

    mImgtbl(output_folder + "/projected", output_folder + "/images.tbl")

    # If selected, perform background matching.

    if background_match:
        
        if verbose:
            print('Calculating overlaps')
        
        mOverlaps(output_folder + "/images.tbl", output_folder + "/diffs.tbl")
        
        if verbose:
            print('Fitting overlap differences')
            
        # If we're properly weighting later this can cause some issues
        # so turn off including the _area files in this case.
        
        mDiffFitExec(
            output_folder + "/projected",
            output_folder + "/diffs.tbl",
            output_folder + "/header.hdr",
            output_folder + "/diffs",
            output_folder + "/fits.tbl",
            noAreas = haveAreas,
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
#             noAreas = haveAreas,
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
        haveAreas=haveAreas,
    )

    # Remove the temp folders we've made along the way

    shutil.rmtree(output_folder + "/projected", ignore_errors=True)
    if background_match:
        shutil.rmtree(output_folder + "/diffs", ignore_errors=True)
        shutil.rmtree(output_folder + "/corrected", ignore_errors=True)
        
def background_median(data,sigma=3,npixels=5,maxiters=20,
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

def make_header(galaxy,width=0.2,height=0.2,
                output_file='hdr.hdr',resolution=1):
    
    """Replacement for mHdr.
    
    mHdr seems to have broken with an IRSA update. This should function
    the same way.
    
    Args:
        galaxy (str): Target to resolve.
        width (float, optional): Output image width in decimal degrees.
            Defaults to 0.2 deg.
        height (float, optional): Output image height in decimal degrees.
            Defaults to 0.2 deg.
        output_file (str, optional): Filename for output header. Defaults
            to hdr.hdr.
        resolution (float, optional): Resolution of output image in 
            arcsec. Defaults to 1.
            
    """
    
    base_url = 'https://irsa.ipac.caltech.edu/cgi-bin/HdrTemplate/nph-hdr?'
    base_url += 'location='+galaxy+'&'
    base_url += 'width='+str(width)+'&'
    base_url += 'height='+str(height)+'&'
    base_url += 'resolution='+str(resolution)
    
    os.system('wget "'+base_url+'"')
    
    # Rename downloaded hdr
    
    hdr_file = glob.glob('nph-hdr*')[0]
    
    os.rename(hdr_file,output_file)

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

def optimize_size(hdu,hdu_out=None):
    
    """Optimize the size of a .fits image.
    
    Calculate the minimum size for an image, updating the HDU at the 
    end. This is useful for images out of mosaicking tools, where there
    might be a number of NaNs artificially inflating the filesize.
    
    Args:
        hdu (str or astropy.io.fits.PrimaryHDU): Data to optimize size of. 
            If this is a string, this is interpreted as a filename for a 
            .fits HDU.
        hdu_out (str, optional): Filename for output trimmed .fits file.
            Defaults to None, which will save nothing.
            
    Returns:
        astropy.io.fits.ImageHDU: The size-optimized HDU.
            
    """
    
    if isinstance(hdu,str):
        hdu = fits.open(hdu)[0]
        data = hdu.data

    wcs = WCS(hdu.header)
        
    x, y = np.meshgrid(range(data.shape[1]),
                       range(data.shape[0]))
    
    idx = np.where(np.isnan(data) == False)
    
    ymin,ymax = np.min(y[idx]),np.max(y[idx])
    xmin,xmax = np.min(x[idx]),np.max(x[idx])
    
    ymid = (ymin+ymax)/2
    xmid = (xmin+xmax)/2
    
    ysize = ymax+1-ymin
    xsize = xmax+1-xmin
    
    cutout = Cutout2D(data,
                      (xmid,ymid),
                      (ysize,xsize),
                      wcs=wcs)
    
    hdu_trimmed = hdu.copy()
    hdu_trimmed.data = cutout.data
    hdu_trimmed.header.update(cutout.wcs.to_header())
    
    if hdu_out is not None:
        hdu_trimmed.writeto(hdu_out,overwrite=True)
    
    return hdu_trimmed

def interp_nans(data,x_stddev=1):
    
    """Interpolate over any NaNs present in a final mosaic.
    
    Uses the Astropy.convolution interpolate_replace_nans to smooth over
    any gaps left in an image. This may be particularly useful for
    WFPC2 images, where there are small gaps between chips.
    
    Args:
        data (numpy.ndarray): Input data to interpolate NaNs over.
        x_stddev (int, optional): Standard deviation of the Gaussian kernel.
            Defaults to 1 (pixel).
            
    Returns:
        numpy.ndarray: The data with NaNs interpolated over
        
    """
    
    kernel = Gaussian2DKernel(x_stddev=x_stddev)
    
    image_interp = interpolate_replace_nans(data,kernel)
    
    return image_interp
