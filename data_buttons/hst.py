# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import fnmatch
import glob
import os
import shutil
import multiprocessing as mp
import logging

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import Observations
from astroquery.gaia import Gaia
from drizzlepac import astrodrizzle, photeq, tweakreg
import stwcs
import crds
import pickle

from . import alignimages

def hst_button(
    galaxies,
    instruments="ACS/WFC",
    prop_ids=None,
    filters=None,
    radius=0.2*u.degree,
    filepath=None,
    download_data=True,
    correct_astrometry=True,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
    log_filename='hst.log',
):
    """Create a HST mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available HST data and mosaics into a final product. It
    will create separate mosaics for each proposal ID, and the file structure
    will look like ``/galaxy/HST/proposal_id/galaxy_instrument_filter_proposal_id.fits``.
    
    N.B. I must confess to not being well-versed with HST data, so if 
    anyone can help improve this please let me know.
    
    This data button uses a number of tools included in the drizzlepac
    Python package. This includes alignimages/tweakreg and astrodrizzle, 
    which correct astrometry and are specifically tailored for the setup 
    of HST data. This means that 1) creating mosaics with this will likely 
    take a long time and 2) you will need a beefy computer (especially with
    regards to hard drive space).
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        instruments (str or list, optional): Instrument to download data 
            for.  Can be any combination of 'ACS/WFC', 'WFC3/IR', 
            'WFC3/UVIS', 'WFPC2/PC', or 'WFPC2/WFC'. If you want all 
            available data for all these instruments, select 'all', but 
            this is not recommended! Defaults to 'ACS/WFC'.
        prop_ids (str or list, optional): Proposal IDs to download data for.
            Defaults to None, which will pull out all proposal IDs for each
            instrument.
        filters (str or list, optional): Filters to download data for.
            The script will look for each filter, for each instrument.
            Defaults to None, which will pull out all applicable filters
            for each instrument, for each proposal ID.
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
        verbose (bool, optional): Can be used to suppress most of the
            output messages produced during the process. Mainly useful
            for debugging. Defaults to False.
        log_filename (str, optional): Will produce a stripped down log
            of what data the code is reducing. By default, will save to
            galaxy/hst.log.
    
    """
    
    if isinstance(galaxies, str):
        galaxies = [galaxies]
        
    if isinstance(instruments,str):
        instruments = [instruments]
        
    if instruments == 'all':
        instruments = ['ACS/WFC',
                       'WFC3/IR','WFC3/UVIS',
                       'WFPC2/PC','WFPC2/WFC']
        
    if isinstance(filters,str):
        filters = [filters]
        
    if isinstance(prop_ids,str):
        prop_ids = [prop_ids]

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
        
        if not os.path.exists(galaxy):
            os.mkdir(galaxy)
            
        if not os.path.exists(galaxy+'/HST'):
            os.mkdir(galaxy+'/HST')
            
        if not verbose:
            
            # Various packages used here put out a lot of messages. Silence info messages.
            
            loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
            for logger in loggers:
                    logger.setLevel(logging.ERROR)
                    
        # Even if verbose is not True, still print out some useful messages to the 
        # console.
        
        logging.basicConfig(filemode='w')
        hst_logger = logging.getLogger('data_buttons')
        handler = logging.FileHandler(galaxy+'/'+log_filename)
        hst_logger.addHandler(handler)
        hst_logger.setLevel(logging.INFO)
        hst_logger.info('Beginning '+galaxy)
        hst_logger.info(' ')
        hst_logger.info(' ')
 
        obs_table = Observations.query_criteria(objectname=galaxy,
                                                radius=radius,
                                                obs_type='all',
                                                obs_collection='HST')
        
        # Ignore any calibration observations.
        obs_table = obs_table[obs_table['intentType'] == 'science']
        
        for instrument in instruments:
    
            # Pixel sizes for final mosaics selected to match the HLA.
            
            pix_size = {'ACS/HRC':0.025,
                        'ACS/SBC':0.03,
                        'ACS/WFC':0.05,
                        'NICMOS/NIC1':0.025,
                        'NICMOS/NIC2':0.05,
                        'NICMOS/NIC3':0.1,
                        'WFC3/IR':0.04,
                        'WFC3/UVIS':0.09,
                        'WFPC2/PC':0.05,
                        'WFPC2/WFC':0.1}[instrument]
                        
            # Filename extension, in order of preference.
            
            suffixes = {'ACS/WFC':['FLC','FLT'],
                        'WFC3/IR':['FLT'],
                        'WFC3/UVIS':['FLC','FLT'],
                        'WFPC2/PC':['C0M'],
                        'WFPC2/WFC':['C0M'],
                        }[instrument]
            
            # The instruments often have / in the name, so account for 
            # this in making folders and files.
            
            hst_logger.info('Beginning '+instrument)
            
            if not os.path.exists(galaxy+'/HST/'+instrument.replace('/','_')):
                os.mkdir(galaxy+'/HST/'+instrument.replace('/','_'))
                
            reset_filters = False
            
            instrument_table = obs_table[obs_table['instrument_name'] == instrument]
            
            reset_prop_ids = False
                
            if not prop_ids:
                prop_ids = list(np.unique(instrument_table['proposal_id']))
                reset_prop_ids = True
                
            hst_logger.info('Available proposal IDs: '+','.join(prop_ids))
            hst_logger.info(' ')
                
            for prop_id in prop_ids:
                
                hst_logger.info('Proposal ID: '+str(prop_id))
                
                prop_table = instrument_table[instrument_table['proposal_id'] == prop_id]
            
                if not filters:
                    filters = list(np.unique(prop_table['filters']))
                    reset_filters = True
                    
                hst_logger.info('Available filters: '+','.join(filters))
            
                for hst_filter in filters:
                    
                    # If we have a highly illegal filter, just skip.
                    # TODO: This needs to be sorted for some fringe
                    # cases, probably.
                     
                    if not hst_filter[0] == 'F':
                        continue
                    
                    hst_logger.info('Filter: '+str(hst_filter))

                    # Pull out available data and download.

                    filter_table = prop_table[prop_table['filters'] == hst_filter]
                    
                    if len(filter_table) == 0:
                        hst_logger.warning('No available data to download. Skipping...')
                        continue
                 
                    data_products_id = Observations.get_product_list(filter_table)
                    
                    for suffix in suffixes:
                    
                        download_table = Observations.filter_products(data_products_id,
                                                                      productSubGroupDescription=suffix,
                                                                      mrp_only=False)
                        
                        if len(download_table) > 0:
                            break
                        
                    filename_ext = suffix.lower()
                    
                    hst_logger.info(instrument+'/'+prop_id+'/'+hst_filter)
                        
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
                     
                    if 1 in steps:
                            
                        # Download files
                        
                        download_mast(download_table,
                                      download_dir="hst_temp/" + galaxy)  
                                    
                        if not os.path.exists(full_filepath+'/raw'):
                            os.mkdir(full_filepath+'/raw')
                        if not os.path.exists(full_filepath+'/outputs'):
                            os.mkdir(full_filepath+'/outputs')    
                                    
                        # Pull out the relevant files, and move to base folder.
                
                        matches = []
                        for root, _, filenames in os.walk("hst_temp/" + galaxy):
                            for filename in fnmatch.filter(
                                filenames, "*_"+filename_ext+".fits"
                            ):
                                matches.append(os.path.join(root, filename))
                
                        for match in matches:
                            
                            filename = match.split('/')
                
                            os.rename(match,full_filepath+'/raw/'+filename[-1])
                            
                        # Clean up any temporary files.
        
                        shutil.rmtree("hst_temp/" + galaxy, ignore_errors=True)
                        
                    hst_files = glob.glob(full_filepath+'/raw/*_'+filename_ext+'.fits')
                        
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
                        
                        # Normalize all files.
                        
                        photeq.photeq(', '.join(hst_files),readonly=False)
                        os.rename('photeq.log','../outputs/photeq.log')
                        
                        if 'WFPC' in instrument:
                            
                            # Using tweakreg, align each frame to GAIA.
                            
                            gaia_table = Gaia.query_object_async(coordinate=galaxy, 
                                      
                                                                 radius=2*radius)
                            ras = gaia_table['ra']
                            decs = gaia_table['dec']
                            
                            source_table = Table([ras,decs])
                            source_table.write('gaia.cat',
                                               format='ascii.fast_commented_header')
                            
                            tweakreg.TweakReg(hst_files,
                                              imagefindcfg={'threshold':5,'conv_width':3},
                                              refcat='gaia.cat',
                                              #expand_refcat=True,
                                              enforce_user_order=False,
                                              shiftfile=True,
                                              outshifts='shifts.txt',
                                              searchrad=10,
                                              minobj=5,
                                              separation=0,
                                              updatehdr=True,
                                              reusename=True,
                                              wcsname='TWEAK',
                                              interactive=False,
                                              fitgeometry='general',
                                              clean=True,
                                              see2dplot=False
                                              )
                            
                            plot_files = glob.glob('*.png')
                            for plot_file in plot_files:
                                os.remove(plot_file)
                                
                            cat_files = glob.glob('*.coo')
                            for cat_file in cat_files:
                                os.remove(cat_file)
                                
                            os.rename('shifts_wcs.fits','../outputs/shifts_wcs.fits')
                            os.rename('tweakreg.log','../outputs/tweakreg.log')
                            os.rename('shifts.txt','../outputs/shifts.txt')
                        
                        elif 'ACS' in instrument or 'WFC3' in instrument:
                             
                            # Correct astrometry using alignimages. First,
                            # correct each frame separately.
                        
                            pool = mp.Pool(mp.cpu_count())
                            
                            suitable_hst_files = pool.map(astrometric_correction,
                                                          hst_files)
                            
                            pool.close()
                            
                            suitable_hst_files = [x for x in suitable_hst_files 
                                                  if x is not None]
                            
                            if len(suitable_hst_files) == 0:
                                hst_logger.warning('Failure with astrometry corrections. Skipping')
                                os.chdir(orig_dir)
                                continue
                            
                            # Now, align every suitable frame simultaneously. 
    
                            output_table = astrometric_correction(suitable_hst_files)
                             
                            with open('../outputs/astrometry.pkl','wb') as table_file:
                                pickle.dump(output_table,table_file)
                            
                        else:
                            
                            raise Exception('Unknown instrument!')
                        
                        os.chdir(orig_dir)
                        
                    os.chdir(full_filepath)
                            
                    if 3 in steps:    
                        
                        os.chdir('raw')
                        
                        if 'WFPC2' in instrument:
                            
                            hst_files = glob.glob('*_c0m.fits')
                            
                            wcskey = 'TWEAK'
                            
                        elif 'ACS' in instrument or 'WFC3' in instrument:
                            
                            with open('../outputs/astrometry.pkl','rb') as table_file:
                             
                                output_table = pickle.load(table_file)
                                    
                            # We only want fits where an acceptable astrometric
                            # solution has been found.
                            
                            suitable_fits = np.where(output_table['fit_qual'] < 5)
#                         (output_table['fit_qual'] >= 1)
                            
                            hst_files = list(output_table[suitable_fits]['imageName'])
                            
                            if len(output_table[suitable_fits]) == 0:
                                hst_logger.warning('Failure with astrometry corrections. Skipping')
                                os.chdir(orig_dir)
                                continue
                            
                            wcskey = ' '
                            
                        else:
                            
                            raise Exception('Unknown instrument!')
        
                        # Perform the mosaicking.
                        
                        #if 'WFPC2' in instrument:
                        #    combine_type = 'median'
                        #elif 'ACS' in instrument:
                        #    combine_type = 'minmed'

                        # Combination types, in order of preference.
                        # Sometimes minmed will fail.

                        combine_types = ['minmed','median']
                        
                        for combine_type in combine_types:
                         
                            try:
                                
                                astrodrizzle.AstroDrizzle(input=hst_files,
                                                          output='../outputs/'+galaxy,
                                                          preserve=False,
                                                          clean=True,
                                                          combine_type=combine_type,
                                                          skymethod='globalmin+match',
                                                          driz_sep_bits='64,32',
                                                          final_scale=pix_size,
                                                          wcskey=wcskey,
                                                          final_rot=0)
                                
                                break
                                
                            except ValueError:
                            
                                pass
                            
                        # Move the AstroDrizzle log.
                        
                        os.rename('astrodrizzle.log',
                                  '../outputs/astrodrizzle.log')
                            
                    # Move back to the original directory.
                        
                    os.chdir(orig_dir)
                    
                    if 4 in steps:
                        
                        # Mosaic outputs, in order of preference.
                        
                        drizzle_exts = ['drc','drz']
                         
                        for drizzle_ext in drizzle_exts:
                         
                            try:
                                
                                convert_to_jy(full_filepath+'/outputs/'+galaxy+'_'+drizzle_ext+'_sci.fits',
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
                                
                                break
                                
                            except IOError:
                                
                                pass
                            
                if reset_filters:
                    filters = None
                    
                hst_logger.info(' ')
        
            if reset_prop_ids:
                prop_ids = None
                
            hst_logger.info(' ')
                    
def download_mast(obs_table,**kwargs):
    
    """Download data from the MAST archive.
    
    This function gives you an easy way to download data from the MAST
    archive, if you already have a table of observations (from Astroquery)
    to hand. By default, this will download everything into your current
    working directory.
    
    Args:
        obs_table (astropy.table.Table): The result of an astroquery.mast
            query.
    
    """
                  
    Observations.download_products(obs_table,**kwargs
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
                                             print_fit_parameters=False,
                                             update_hdr_wcs=True,
                                             )
    
    os.rename(runfile,
              '../outputs/'+runfile)
                            
    if single_file:
        
        # Filter the table down to only one where sources have been detected
        # and the fit is not totally broken. Else alignimages on the whole
        # set can cause issues.
        
        if output_table['fit_qual'] < 5 and output_table['foundSources'] > 0:
#             1 <= 
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
        