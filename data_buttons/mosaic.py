# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os
import shutil

from MontagePy.archive import *
from MontagePy.main import *


def mosaic(
        input_folder,
        header=None,
        output_folder='mosaic',
        background_match=True,
):


    """Mosaic together a folder full of .fits files.
    
    Takes a folder and runs through the number of steps that MontagePy
    requires to combine all the .fits files within that folder. Also,
    optionally performs background modelling and matching (which most
    users will want to set to true). Note that this is not a background
    subtraction!
    
    Args:
        input_folder (str): Folder of raw files to mosaic.
        header (str, optional): Output from mHdr. If not specified, 
            will mosaic all of the images together not matter the 
            overlap.
        output_folder (str, optional): Working folder for mosaicking.
            Defaults to 'mosaic'.
        background_match (bool, optional): Whether to perform background
            matching steps while mosaicking. Defaults to True.
            
    Todo:
        * Create verbose print statements, for debugging.
            
    """

    # Create the working folders we'll need

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    os.mkdir(output_folder + '/projected')
    if background_match:
        os.mkdir(output_folder + '/diffs')
        os.mkdir(output_folder + '/corrected')

    # Make an optimum header for these images

    _ = mImgtbl(
        input_folder,
        output_folder + '/images.tbl',
        )
    
    if header is None:
        _ = mMakeHdr(
            output_folder + '/images.tbl',
            output_folder + '/header.hdr',
            )
    else:
        shutil.copy(header,
                    output_folder + '/header.hdr')

    # Project the original images to this header

    _ = mProjExec(
        input_folder,
        output_folder + '/images.tbl',
        output_folder + '/header.hdr',
        projdir=output_folder + '/projected',
        quickMode=True)

    _ = mImgtbl(
        output_folder + '/projected',
        output_folder + '/images.tbl',
        )

    # If selected, perform background matching.

    if background_match:
        _ = mOverlaps(
            output_folder + '/images.tbl',
            output_folder + '/diffs.tbl',
            )
        _ = mDiffFitExec(
            output_folder + '/projected',
            output_folder + '/diffs.tbl',
            output_folder + '/header.hdr',
            output_folder + '/diffs',
            output_folder + '/fits.tbl',
            )
        _ = mBgModel(
            output_folder + '/images.tbl',
            output_folder + '/fits.tbl',
            output_folder + '/corrections.tbl',
            )
        _ = mBgExec(
            output_folder + '/projected',
            output_folder + '/images.tbl',
            output_folder + '/corrections.tbl',
            output_folder + '/corrected',
            )
        _ = mImgtbl(
            output_folder + '/corrected',
            output_folder + '/images.tbl')

    # Finally, coadd the images

    if not background_match:
        folder = output_folder + '/projected'
    else:
        folder = output_folder + '/corrected'

    _ = mAdd(
        folder,
        output_folder + '/images.tbl',
        output_folder + '/header.hdr',
        output_folder + '/mosaic.fits',
        )

    # Remove the temp folders we've made along the way

    shutil.rmtree(output_folder + '/projected',
                  ignore_errors=True,
                  )
    if background_match:
        shutil.rmtree(output_folder + '/diffs',
                      ignore_errors=True,
                      )
        shutil.rmtree(output_folder + '/corrected',
                      ignore_errors=True,
                      )
