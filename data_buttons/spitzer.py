# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import fnmatch
import glob
import os
import shutil

import astropy.units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astroquery.ned import Ned
from astroquery.simbad import Simbad
from MontagePy.main import mProject, mHdr
import numpy as np
import wget

from . import tools


def spitzer_button(
    galaxies,
    filters="all",
    radius=None,
    filepath=None,
    download_data=True,
    create_mosaic=True,
    jy_conversion=True,
    verbose=False,
):
    """Create an Spitzer mosaic, given a galaxy name.
    
    Using a galaxy name and radius, queries around that object, 
    downloads available Spitzer data from the SHA and mosaics.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
        filters (str or list, optional): Any combination of 'IRAC1', 'IRAC2',
            'IRAC3', 'IRAC4', 'MIPS24', 'MIPS70', or 'MIPS160'. If you 
            want everything, select 'all'. Defaults to 'all'.
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
    
    base_url = 'https://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/DataService'

    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if filters == "all":
        filters = ['IRAC1','IRAC2','IRAC3','IRAC4',
                   'MIPS24','MIPS70','MIPS160']

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
        
    if not os.path.exists('Spitzer_temp'):
        os.mkdir('Spitzer_temp')

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
            
        if not os.path.exists(galaxy + "/Spitzer/"):
            os.mkdir(galaxy + "/Spitzer/")
            
        # Download VOTable covering the user-specified area. Resolve
        # position of galaxy
        
        simbad_query = Simbad.query_object(galaxy)
        
        ra = simbad_query[0]['RA']
        dec = simbad_query[0]['DEC']
        
        # Convert to decimal degrees
        
        c = SkyCoord(ra=ra,dec=dec,unit=(u.hourangle,u.deg))
        
        diam_deg = str(2*radius.to(u.degree).value)
        
        query_url = base_url+'?RA='+str(c.ra.value)+'&DEC='+str(c.dec.value)+'&SIZE='+str(diam_deg)
        query_url += '&VERB=3&DATASET=ivo%3A%2F%2Firsa.csv%2Fspitzer.level2'
        
        if os.path.exists(galaxy+'/Spitzer/table.xml'):
            os.remove(galaxy+'/Spitzer/table.xml')

        wget.download(query_url,
                      out=galaxy+'/Spitzer/table.xml')
        
        votable = ascii.read(galaxy+'/Spitzer/table.xml')

        for spitzer_filter in filters:
            
            # Name of wavelengths for the table
            
            table_wl = {'IRAC1':'IRAC 3.6um',
                        'IRAC2':'IRAC 4.5um',
                        'IRAC3':'IRAC 5.8um',
                        'IRAC4':'IRAC 8.0um',
                        'MIPS24':'MIPS 24um',
                        'MIPS70':'MIPS 70um',
                        'MIPS160':'MIPS 160um'}[spitzer_filter]
            
            # Output mosaic pixel scales in arcsec
            
            pixel_scale = {'IRAC1':0.6,
                           'IRAC2':0.6,
                           'IRAC3':0.6,
                           'IRAC4':0.6,
                           'MIPS24':2.45,
                           'MIPS70':4,
                           'MIPS160':8}[spitzer_filter]
            
            if verbose:
                print('Beginning '+spitzer_filter)

            if not os.path.exists(galaxy + "/Spitzer/" + spitzer_filter):
                os.mkdir(galaxy + "/Spitzer/" + spitzer_filter)
                
            if not os.path.exists(galaxy + "/Spitzer/" + spitzer_filter+'/raw'):
                os.mkdir(galaxy + "/Spitzer/" + spitzer_filter+'/raw')
                
            if not os.path.exists(galaxy + "/Spitzer/" + spitzer_filter+'/data'):
                os.mkdir(galaxy + "/Spitzer/" + spitzer_filter+'/data')
            
            if not os.path.exists(galaxy + "/Spitzer/" + spitzer_filter+'/weight'):
                os.mkdir(galaxy + "/Spitzer/" + spitzer_filter+'/weight')
            
            if not os.path.exists(galaxy + "/Spitzer/" + spitzer_filter+'/outputs'):
                os.mkdir(galaxy + "/Spitzer/" + spitzer_filter+'/outputs')
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
                    
                if not os.path.exists('Spitzer_temp/'+galaxy):
                    os.mkdir('Spitzer_temp/'+galaxy)
    
                # Match up the wavelength column of the table to the 
                # correct name.
                
                for row in votable:
                               
                    if row['wavelength'] == table_wl:
                                   
                        # Download files along with ancillary stuff.
                                   
                        dl_file = row['accessWithAnc1Url']
                        
                        # Sometimes we get a NONE so in that case, skip.
                        
                        if dl_file == 'NONE':
                            continue
                                   
                        wget.download(
                            dl_file,
                            out='Spitzer_temp/'+galaxy+'/')
                                   
                # Unzip all these files
                           
                spitzer_files = glob.glob('Spitzer_temp/'+galaxy+'/*.zip')
                           
                for spitzer_file in spitzer_files:
                    os.system('unzip '+spitzer_file+' -d '+'Spitzer_temp/'+galaxy)
                               
                # Move any maic and munc files to raw
                           
                extensions = ['maic','munc']
                           
                for extension in extensions:
                            
                    for root, _, filenames in os.walk("Spitzer_temp/" + galaxy):
                        for filename in fnmatch.filter(
                            filenames, "*"+extension+".fits"
                        ):
                            match = os.path.join(root, filename)
                                        
                            filename = match.split('/')
                                        
                            os.rename(
                                    match,
                                    galaxy
                                    + "/Spitzer/"
                                    + spitzer_filter
                                    + "/raw/"
                                    + filename[-1],
                                )
                            
                shutil.rmtree('Spitzer_temp/'+galaxy,ignore_errors=True)
                 
                spitzer_files = glob.glob(galaxy + "/Spitzer/" + spitzer_filter+'/raw/*maic.fits')
                            
                # Reproject, taking into account the uncertainties to weight
                # the maps initially
                 
                if verbose:
                    print('Performing inital round of weighted reprojection')
                 
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/Spitzer/"+spitzer_filter+"/outputs/header.hdr",
                    resolution=pixel_scale,
                    )
                 
                for spitzer_file in spitzer_files:
                     
                    # Do a little fiddling to find the name of the 
                    # uncertainty file.
                     
                    unc_name = spitzer_file.split('/')[-1]
                    unc_name = '_'.join(unc_name.split('_')[2:5])
                    unc_file = glob.glob(galaxy + "/Spitzer/" + spitzer_filter+'/raw/*'+unc_name+'*munc.fits')
                    unc_file = unc_file[0]
                     
                    # Reproject each image to final pixel scale, weighting 
                    # by inverse uncertainty.
                     
                    hdu_exp = fits.open(unc_file)[0]
                    hdu_exp.data = hdu_exp.data**-1
                     
                    hdu_exp.writeto(unc_file.replace('raw','weight'),
                                    overwrite=True)
                     
                    mProject(spitzer_file,
                             spitzer_file.replace('raw','data'),
                             galaxy + "/Spitzer/"+spitzer_filter+"/outputs/header.hdr",
                             weight_file=unc_file.replace('raw','weight'))
                
            if 2 in steps:
    
                # Mosaic all these files together.
    
                if verbose:
                    print("Beginning mosaics")
                    
                mHdr(
                    galaxy,
                    2 * radius.value,
                    2 * radius.value,
                    galaxy + "/Spitzer/"+spitzer_filter+"/outputs/header.hdr",
                    resolution=pixel_scale,
                    )
 
                tools.mosaic(
                    galaxy + "/Spitzer/" + spitzer_filter+'/data',
                    header=galaxy+"/Spitzer/"+spitzer_filter+"/outputs/header.hdr",
                    verbose=verbose,
                    reproject=False,
                    haveAreas=True,
                    )
                   
                os.rename(
                    "mosaic/mosaic.fits", 
                    galaxy + "/Spitzer/"+spitzer_filter+"/outputs/"+galaxy+'.fits'
                    )
 
                shutil.rmtree("mosaic/", ignore_errors=True)
                
            if 3 in steps:
                
                if verbose:
                    print('Converting to Jy')
            
                # Convert to Jy.
                
                convert_to_jy(galaxy + "/Spitzer/" + spitzer_filter + "/outputs/"+galaxy+".fits",
                              hdu_out=galaxy + "/Spitzer/"+galaxy+"_" + spitzer_filter+".fits")
                
        if original_radius is None:
            radius = None
        else:
            radius = original_radius.copy()
            
def convert_to_jy(hdu_in,hdu_out=None):
    
    """Convert from MJy/sr to Jy/px.
    
    Spitzer maps are in MJy/sr, which can be converted through to units
    of Jy/px.
    
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of Spitzer
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
