# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import glob
import os
import shutil

import astropy.units as u
from astropy.io import fits
from astroquery.ned import Ned
from astroquery.esasky import ESASky
from MontagePy.main import mProject, mHdr
import numpy as np

from . import tools


def herschel_button(
    galaxies,
    filters="all",
    radius=None,
    filepath=None,
    download_data=True,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
):
    """Create an Herschel mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available Herschel data from the HSA and mosaics.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
        filters (str or list, optional): Any combination of 'PACS70',
            'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500'. If 
            you want everything, select 'all'. Defaults to 'all'.
        radius (astropy.units.Quantity, optional): Radius around the 
            galaxy to search for observations. Defaults to None, where
            it will query Ned to get size.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        download_data (bool, optional): If True, will download data using 
            Astroquery. Defaults to True.
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

    if filters == "all":
        filters = ['PACS70','PACS100','PACS160',
                   'SPIRE250','SPIRE350','SPIRE500']

    if isinstance(filters, str):
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
        
    if not os.path.exists('Herschel_temp'):
        os.mkdir('Herschel_temp')
        
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
            
        if not os.path.exists(galaxy + "/Herschel/"):
            os.mkdir(galaxy + "/Herschel/")
            
        for herschel_filter in filters:
            
            pixel_scale = {'PACS70':2,
                           'PACS100':3,
                           'PACS160':4,
                           'SPIRE250':6,
                           'SPIRE350':8,
                           'SPIRE500':12}[herschel_filter]
            
            if not os.path.exists(galaxy + "/Herschel/"+herschel_filter):
                os.mkdir(galaxy + "/Herschel/"+herschel_filter)
            if not os.path.exists(galaxy + "/Herschel/"+herschel_filter+"/raw"):
                os.mkdir(galaxy + "/Herschel/"+herschel_filter+"/raw")
            if not os.path.exists(galaxy + "/Herschel/"+herschel_filter+"/data"):
                os.mkdir(galaxy + "/Herschel/"+herschel_filter+"/data")
            if not os.path.exists(galaxy + "/Herschel/"+herschel_filter+"/reprojected"):
                os.mkdir(galaxy + "/Herschel/"+herschel_filter+"/reprojected")
            if not os.path.exists(galaxy + "/Herschel/"+herschel_filter+"/weight"):
                os.mkdir(galaxy + "/Herschel/"+herschel_filter+"/weight")
            if not os.path.exists(galaxy + "/Herschel/"+herschel_filter+"/outputs"):
                os.mkdir(galaxy + "/Herschel/"+herschel_filter+"/outputs")
            
            if 1 in steps:
                
                if verbose:
                    print('Downloading available data')
                 
                # Download available Herschel images for the object.
              
                images = ESASky.get_images(galaxy, radius=0.4*u.degree,
                                           missions=['Herschel'],
                                           download_dir='Herschel_temp/'+galaxy)
                  
                images = images['HERSCHEL']
                  
                # Pull out available data for each waveband and save.
                  
                herschel_key = {'PACS70':'70',
                                'PACS100':'100',
                                'PACS160':'160',
                                'SPIRE250':'250',
                                'SPIRE350':'350',
                                'SPIRE500':'500'}[herschel_filter]
                                  
                i = 0
          
                for image in images:
                    for key in image.keys():
                        if key == herschel_key:
                              
                            image[key].writeto(galaxy + "/Herschel/"+herschel_filter+"/raw/"+herschel_filter+"_"+str(i)+".fits",
                                               overwrite=True)
                              
                            print(image[key][0].header['INSTMODE'])
                              
                            i += 1
                          
                # We now want to pull out the data and the coverage map for 
                # each image.
                  
                if verbose:
                    print('Beginning initial weighted reprojection')
                  
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/Herschel/"+herschel_filter+"/outputs/header.hdr",
                    resolution=pixel_scale,
                    )
                   
                herschel_files = glob.glob(galaxy + "/Herschel/"+herschel_filter+"/raw/*.fits")
                   
                for herschel_file in herschel_files:
                       
                    hdu_data = fits.open(herschel_file)[1]
                    hdu_data.data[hdu_data.data == 0] = np.nan
                    
                    orig_pix_scale = np.abs(hdu_data.header['CDELT1'])
                                           
                    # Coverage is proportional to the exposure time so 
                    # use sqrt of that as the weight.
                       
                    hdu_weight = fits.open(herschel_file)[2]
                    hdu_weight.data[hdu_weight.data == 0] = np.nan
                    hdu_weight.data = hdu_weight.data**0.5
                       
                    fits.writeto(herschel_file.replace('/raw/','/data/'),
                                 hdu_data.data,hdu_data.header,
                                 overwrite=True)
                    fits.writeto(herschel_file.replace('/raw/','/weight/'),
                                 hdu_weight.data,hdu_weight.header,
                                 overwrite=True)
                      
                    # Perform an initial reprojection, weighting by the
                    # weight map
                      
                    mProject(herschel_file.replace('raw','data'),
                             herschel_file.replace('raw','reprojected'),
                             galaxy + "/Herschel/"+herschel_filter+"/outputs/header.hdr",
                             weight_file=herschel_file.replace('raw','weight'),
                             )
                    
                    try:
                        
                        hdu_reproj = fits.open(herschel_file.replace('raw','reprojected'))[0]
                        
                        # If PACS, account for pixel scale
                         
                        if 'PACS' in herschel_filter:
                              
                            new_pix_scale = np.abs(hdu_reproj.header['CDELT1'])
                              
                            hdu_reproj.data *= new_pix_scale**2/orig_pix_scale**2
                        
                        hdu_reproj.writeto(herschel_file.replace('raw','reprojected'),
                                           overwrite=True)
                        
                    except FileNotFoundError:
                        
                        pass
                    
            if 2 in steps:
                
                # Mosaic everything together.
                
                if verbose:
                    print("Beginning mosaics")
                    
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/Herschel/"+herschel_filter+"/outputs/header.hdr",
                    resolution=pixel_scale,
                    )
 
                tools.mosaic(
                    galaxy + "/Herschel/" + herschel_filter+'/reprojected',
                    header=galaxy+"/Herschel/"+herschel_filter+"/outputs/header.hdr",
                    verbose=verbose,
                    reproject=False,
                    haveAreas=True,
                    )
                   
                os.rename(
                    "mosaic/mosaic.fits", 
                    galaxy + "/Herschel/"+herschel_filter+"/outputs/"+galaxy+'.fits'
                    )
 
                shutil.rmtree("mosaic/", ignore_errors=True)
                
            if 3 in steps:
                
                # Convert to Jy.
                
                if verbose:
                    print('Converting to Jy')
                    
                # For PACS, we already have this in Jy so
                
                if 'PACS' in herschel_filter:
                
                    shutil.copy(galaxy + "/Herschel/"+herschel_filter+"/outputs/"+galaxy+'.fits',
                                galaxy + "/Herschel/"+galaxy+'_'+herschel_filter+'.fits')
                
                # For SPIRE, this is in MJy/sr
                
                if 'SPIRE' in herschel_filter:
                    
                    convert_to_jy(galaxy + "/Herschel/"+herschel_filter+"/outputs/"+galaxy+'.fits',
                                  galaxy + "/Herschel/"+galaxy+'_'+herschel_filter+'.fits')
                        
        shutil.rmtree('Herschel_temp/'+galaxy,
                      ignore_errors=True)
        if original_radius is None:
            radius = None
        else:
            radius = original_radius.copy()
            
def convert_to_jy(hdu_in,hdu_out=None):
    
    """Convert from MJy/sr to Jy/px.
    
    SPIRE maps are in MJy/sr, which can be converted through to units
    of Jy/px.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of SPIRE
            .fits file, or an Astropy PrimaryHDU instance (i.e. the result 
            of ``fits.open(file)[0]``).
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
        
    # Convert from MJy/sr to Jy/px
    
    data *= 1e6 #MJy to Jy
    
    steradian_to_arcsec = 4.25e10
    pix_arcsec = np.abs(header['CDELT1'])*3600
    
    data *= pix_arcsec**2/steradian_to_arcsec
     
    header['BUNIT'] = 'Jy/pix'
    
    if hdu_out is not None:
        fits.writeto(hdu_out,
                     data,header,
                     overwrite=True)
        
    return fits.PrimaryHDU(data=data,header=header)
    