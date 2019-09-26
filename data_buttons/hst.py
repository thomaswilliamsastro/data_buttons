# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import fnmatch
import glob
import os
import shutil
import sys
import multiprocessing as mp

import astropy.units as u
import numpy as np
from astropy.io import fits
from astroquery.mast import Observations
from drizzlepac import astrodrizzle, photeq
import stwcs
import crds
import pickle

from . import alignimages

def hst_button(
    galaxies,
    instruments="ACS/WFC",
    filters="all",
    radius=0.2*u.degree,
    filepath=None,
    download_data=True,
    correct_astrometry=True,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
):
    """Create a HST mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available HST data and mosaics into a final product. It
    will create separate mosaics for each proposal ID, and the file structure
    will look like ``/galaxy/HST/proposal_id/galaxy_instrument_filter_proposal_id.fits``.
    
    N.B. This currently has been tested only with ACS/WFC data. I must
    also confess to not being well-versed with HST data, so if anyone
    can help improve this please let me know.
    
    This data button uses a number of tools included in the drizzlepac
    Python package. This includes alignimages and astrodrizzle, which
    correct astrometry and are specifically tailored for the setup of HST
    data. This means that 1) creating mosaics with this will likely take
    a long time and 2) you will need a beefy computer (especially with
    regards to hard drive space).
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        instruments (str or list): Instrument to download data for. 
            Currently only works for ACS/WFC. Defaults to 'ACS/WFC'.
        filters (str or list, optional): Filters to download data for.
            The script will look for each filter, for each instrument.
            Defaults to 'all'.
        radius (astropy.units.Quantity, optional): Radius around the 
            galaxy to search for observations. Defaults to 0.2*u.degrees.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        download_data (bool, optional): If True, will download data from 
            MAST. Defaults to True.
        correct_astrometry (bool, optional): If True, will perform astrometric
            corrections to the downloaded data using alignimages. Defaults
            to True.
        create_mosaic (bool, optional): Switching this to True will 
            mosaic the data using astrodrizzle as appropriate. Defaults 
            to True.
        jy_conversion (bool, optional): Convert the mosaicked file from
            raw units to Jy/pix. Defaults to True.
        verbose (bool, optional): Print out messages during the process.
            Useful mainly for debugging purposes or large images. 
            Defaults to False.
            
    Todo:
        * Extend to other HST instruments.
    
    """
    
    if isinstance(galaxies, str):
        galaxies = [galaxies]
        
    if instruments != 'ACS/WFC':
        raise Exception('This tool currently only works for ACS/WFC.')
        
    if isinstance(instruments,str):
        instruments = [instruments]
        
    if filters == 'all':
        filters = ['F220W','F250W','F330W',
                   'F344N','F435W','F475W',
                   'F502N','F555W','F550M',
                   'F606W','F625W','F658N',
                   'F660N','F775W','F814W',
                   'F892N','F850LP']

    if isinstance(filters,str):
        filters = [filters]

    if filepath is not None:
        os.chdir(filepath)
        orig_dir = filepath
    else:
        orig_dir = os.getcwd()
        
    steps = []
    
    if download_data:
        steps.append(1)
    if correct_astrometry:
        steps.append(2)
    if create_mosaic:
        steps.append(3)
    if jy_conversion:
        steps.append(4)
        
    # Set up folders for various corrections
    
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH'] = 'reference_files'
    os.environ['iref'] = 'reference_files/references/hst/wfc3/'
    os.environ['jref'] = 'reference_files/references/hst/acs/'
    os.environ['uref'] = 'reference_files/references/hst/wfpc2/'
        
    for galaxy in galaxies:
        
        if verbose:
            print('Beginning '+galaxy)

        if not os.path.exists(galaxy):
            os.mkdir(galaxy)
            
        if not os.path.exists(galaxy+'/HST'):
            os.mkdir(galaxy+'/HST')
 
        obs_table = Observations.query_criteria(objectname=galaxy,
                                                radius=radius,
                                                obs_type='all',
                                                obs_collection='HST')
        
        for instrument in instruments:
            
            # The instruments often have / in the name, so account for 
            # this in making folders and files.
            
            if verbose:
                print('Beginning '+instrument)
                
            
            if not os.path.exists(galaxy+'/HST/'+instrument.replace('/','_')):
                os.mkdir(galaxy+'/HST/'+instrument.replace('/','_'))
            
            for hst_filter in filters:
                    
                prop_ids = np.unique(obs_table['proposal_id'])
                
                for prop_id in prop_ids:
                    
                    # Pull out available data and download
    
                    query_results = np.where((obs_table['instrument_name'] == instrument) & 
                                             (obs_table['filters'] == hst_filter) & 
                                             (obs_table['proposal_id'] == prop_id) &
                                             (obs_table['intentType'] == 'science'))
                 
                    # If there isn't any HST coverage, just skip
         
                    if len(query_results[0]) == 0:
                        continue
                    
                    if verbose:
                        print('Beginning '+instrument+' '+hst_filter+' Proposal ID '+prop_id)
                        
                    if not os.path.exists(galaxy+
                                          '/HST/'+
                                          instrument.replace('/','_')+
                                          '/'+
                                          hst_filter):
                        os.mkdir(galaxy+
                                 '/HST/'+
                                 instrument.replace('/','_')+
                                 '/'+
                                 hst_filter) 
                        
                    if not os.path.exists(galaxy+'/HST/'+prop_id):
                        os.mkdir(galaxy+'/HST/'+prop_id) 
                        
                    full_filepath =  (galaxy+
                                          '/HST/'+
                                          instrument.replace('/','_')+
                                          '/'+
                                          hst_filter+
                                          '/'
                                          +prop_id)
                    
                    if not os.path.exists(full_filepath):
                        os.mkdir(full_filepath)
                     
                    # We only want to print out download messages if
                    # verbose is True, so redirect otherwise.
                     
                    if 1 in steps:
                         
                        if not verbose:
                            sys.stdout = open(os.devnull,'w')
                            
                        download_mast(obs_table[query_results],
                                      download_dir="hst_temp/" + galaxy,
                                      productSubGroupDescription=['FLC','FLT'])
                                     
                        # And set back to the original for printing.
                                 
                        if not verbose:
                                     
                            sys.stdout = sys.__stdout__    
                                    
                        if not os.path.exists(full_filepath+'/raw'):
                            os.mkdir(full_filepath+'/raw')
                        if not os.path.exists(full_filepath+'/outputs'):
                            os.mkdir(full_filepath+'/outputs')    
                                    
                        # Pull out the relevant files, and move to base folder.
                        # TODO: This will be different for different instruments.
                
                        matches = []
                        for root, _, filenames in os.walk("hst_temp/" + galaxy):
                            for filename in fnmatch.filter(
                                filenames, "*_fl?.fits"
                            ):
                                matches.append(os.path.join(root, filename))
                
                        for match in matches:
                            
                            filename = match.split('/')
                
                            os.rename(match,full_filepath+'/raw/'+filename[-1])
                            
                        # Clean up any temporary files.
        
                        shutil.rmtree("hst_temp/" + galaxy, ignore_errors=True)
                        
                    # Filename extensions, in order of preference
                    
                    filename_exts = {'ACS/WFC':['flc','flt']}[instrument]
                     
                    for filename_ext in filename_exts:
                        
                        hst_files = glob.glob(full_filepath+'/raw/*_'+filename_ext+'.fits')
                        
                        if len(hst_files) > 0:
                            
                            # Remove any extraneous files to keep space used to a minimum.
                            
                            filename_exts.remove(filename_ext)
                            
                            for supercede_ext in filename_exts:
                                for filename in glob.glob(full_filepath+'/raw/*_'+supercede_ext+'.fits'):
                                    os.remove(filename)
                            
                            break
                        
                    if 2 in steps:
                            
                        # First, update the WCS information in case it's 
                        # required.
                           
                        crds.assign_bestrefs(hst_files,
                                             sync_references=True)
            
                        for hst_file in hst_files:
                         
                            stwcs.updatewcs.updatewcs(hst_file,
                                                      use_db=False)
                            
                        os.chdir(full_filepath+'/raw')
                            
                        hst_files = glob.glob('*_'+filename_ext+'.fits')
                             
                        # Correct astrometry using alignimages. First,
                        # correct each frame separately.
                        
                        pool = mp.Pool(mp.cpu_count())
                        
                        suitable_hst_files = pool.map(astrometric_correction,
                                                      hst_files)
                        
                        pool.close()
                        
                        try:
                            suitable_hst_files.remove(None)
                        except ValueError:
                            pass
                        
                        if len(suitable_hst_files) == 0:
                            print('Failure with astrometry corrections. Skipping')
                            os.chdir(orig_dir)
                            continue
                        
                        # Now, align every suitable frame simultaneously. 

                        output_table = astrometric_correction(suitable_hst_files)
                         
                        with open('../outputs/astrometry.pkl','wb') as table_file:
                            pickle.dump(output_table,table_file)
                           
                    else:
                        
                        os.chdir(full_filepath+'/raw')
                         
                        with open('../outputs/astrometry.pkl','rb') as table_file:
                             
                            output_table = pickle.load(table_file)
                            
                    if 3 in steps:    
                                
                        # We only want fits where an acceptable astrometric
                        # solution has been found.
                        
                        suitable_fits = np.where(output_table['fit_qual'] < 5)
#                         (output_table['fit_qual'] >= 1)
                        
                        hst_files = list(output_table[suitable_fits]['imageName'])
                        
                        if len(output_table[suitable_fits]) == 0:
                            print('Failure with astrometry corrections. Skipping')
                            os.chdir(orig_dir)
                            continue
                        
                        # Normalize all suitable files.
                        
                        photeq.photeq(', '.join(hst_files),readonly=False)
                        
                        # And perform the mosaicking.
                        
                        pix_size = {'ACS/WFC':0.05}[instrument]
                         
                        try:
                            
                            astrodrizzle.AstroDrizzle(input=hst_files,
                                                      output='../outputs/'+galaxy,
                                                      preserve=False,
                                                      clean=True,
                                                      skymethod='globalmin+match',
                                                      driz_sep_bits='64,32',
                                                      final_scale=pix_size,
                                                      final_rot=0)
                        except ValueError:
                            
                            # Occasionally, the median backgrounding will
                            # fail. In this case, change to 'median' 
                            # rather than 'minmed' and try again.
                            
                            astrodrizzle.AstroDrizzle(input=hst_files,
                                                  output='../outputs/'+galaxy,
                                                  preserve=False,
                                                  clean=True,
                                                  combine_type='median',
                                                  driz_sep_bits='64,32',
                                                  final_scale=pix_size,
                                                  final_rot=0)
                            
                        # Move the AstroDrizzle log.
                        
                        os.rename('astrodrizzle.log',
                                  '../outputs/astrodrizzle.log')
                            
                        # Move back to the original directory.
                            
                        os.chdir(orig_dir)
                    
                    if 4 in steps:
                         
                        try:
                            convert_to_jy(full_filepath+'/outputs/'+galaxy+'_drc_sci.fits',
                                          galaxy+
                                          '/HST/'
                                          +prop_id
                                          +'/'
                                          +galaxy
                                          +'_'
                                          +instrument.replace('/','_')
                                          +'_'
                                          +hst_filter
                                          +'_'
                                          +prop_id
                                          +'.fits')
                        except:
                            convert_to_jy(full_filepath+'/outputs/'+galaxy+'_drz_sci.fits',
                                          galaxy+
                                          '/HST/'
                                          +prop_id
                                          +'/'
                                          +galaxy
                                          +'_'
                                          +instrument.replace('/','_')
                                          +'_'
                                          +hst_filter
                                          +'_'
                                          +prop_id
                                          +'.fits')
                    
def download_mast(obs_table,
                  **kwargs):
    
    """Download data from the MAST archive.
    
    This function gives you an easy way to download data from the MAST
    archive, if you already have a table of observations (from Astroquery)
    to hand. By default, this will download everything that satisfies the
    information in the observation table, but **kwargs can pare down that
    data to cut down on space.
    
    Args:
        obs_table (astropy.table.Table): The result of an astroquery.mast
            query.
    
    """
            
    dataProductsByID = Observations.get_product_list(obs_table)
                  
    Observations.download_products(
        dataProductsByID, mrp_only=False, 
        **kwargs
    )
        
def astrometric_correction(filename):
    
    """Perform astrometric corrections on HST data.
    
    Using the astrodrizzle alignimages package, will perform astrometric
    corrections on either a single file or a set of files. If a single
    file, then it will only provide absolute astrometric corrections
    based on Gaia observations. If multiple files, it will also correct
    for relative corrections between observations.
    
    Args:
        filename (str or list): The filenames for the data to correct.
            Can either be a single string or a list of filenames. 
            IMPORTANT: These filenames should be just the filename with
            no directory structure beforehand. Else alignimages gets 
            confused!
            
    Returns:
        str or astropy.table.Table: 
            The filename for which an astrometric solution has been found 
            (if input is a string), or information regarding all of the 
            astrometric fits (if input is a list).
    
    """
    
    single_file = False
    
    if isinstance(filename,str):
        single_file = True
        runfile = 'alignimages_'+filename+'.log'
        filename_list = [filename]
    else:
        filename_list = filename.copy()
        runfile = 'alignimages.log'
    
    output_table = alignimages.perform_align(filename_list,
                                             runfile=runfile,
                                             update_hdr_wcs=True,
                                             )
    
    os.rename(runfile,
              '../outputs/'+runfile)
                            
    if single_file:
        if output_table['fit_qual'] < 5:
            return filename
        
    else:                           
        return output_table                 
                    
def convert_to_jy(hdu_in,hdu_out=None):
    
    """Convert from HST electrons/s to Jy/pixel.
    
    HST maps are provided in units of electrons/s. These can be converted
    to more helpful units using the keywords in the .fits header, specifically
    PHOTFLAM (the sensitivity in ergs/cm^2/Ang/electron), and the pivot
    wavelength PHOTPLAM.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of HST 
            .fits file, or an Astropy PrimaryHDU instance (i.e. the 
            result of ``fits.open(file)[0]``).
        hdu_out (str, optional). If not None, will save the hdu to the
            string provided. Defaults to None.
        
    Returns:
        astropy.io.fits.PrimaryHDU: The HDU in units of Jy/pix.
    
    """
    
    if isinstance(hdu_in,str):
        hdu = fits.open(hdu_in)[0]
    else:
        hdu = hdu_in.copy()
    
    data = hdu.data.copy()
    header = hdu.header.copy()
    
    data *= header['PHOTFLAM']
    data *= 3.33564095e4 * header['PHOTPLAM']**2
    
    header['BUNIT'] = 'Jy/pix'
    
    if hdu_out is not None:
        fits.writeto(hdu_out,
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
        