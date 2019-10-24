# Ensure python3 compatibility
from __future__ import absolute_import, print_function, division

import os

import astropy.units as u
from astropy.io import fits
import numpy as np
from astroquery.skyview import SkyView
from astroquery.ned import Ned


def planck_button(
    galaxies,
    frequencies="all",
    radius=None,
    filepath=None,
    download_data=True,
    jy_conversion=True,
    verbose=False,
):
    """Obtain Planck cutout, given a galaxy name.
    
    Using a galaxy name and radius, uses SkyView to download a Planck
    cutout.
    
    Args:
        galaxies (str or list): Names of galaxies to create mosaics for.
            Resolved by NED.
        frequencies (str or list, optional): Any combination of '030',
            '044', '070', '100', '143', '217', '353', '545', '857'. If you 
            want everything, select 'all'. Defaults to 'all'.
        radius (astropy.units.Quantity, optional): Radius around the 
            galaxy to search for observations. Defaults to None, where
            it will query Ned to get size.
        filepath (str, optional): Path to save the working and output
            files to. If not specified, saves to current working 
            directory.
        download_data (bool, optional): If True, will download data from
            SkyView. Defaults to True.
        jy_conversion (bool, optional): Convert the downloaded file from
            raw units to Jy/pix. Defaults to True.
        verbose (bool, optional): Print out messages during the process.
            Useful mainly for debugging purposes or large images. 
            Defaults to False.
    
    """

    if isinstance(galaxies, str):
        galaxies = [galaxies]

    if frequencies == "all":
        frequencies = ['030','044','070','100','143','217','353','545','857']

    if isinstance(frequencies, str):
        frequencies = [frequencies]

    if filepath is not None:
        os.chdir(filepath)
        
    if radius is not None:
        original_radius = radius.copy()
    else:
        original_radius = None
        
    steps = []
    
    if download_data:
        steps.append(1)
    if jy_conversion:
        steps.append(2)

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

        for planck_freq in frequencies:
            
            # Generate cutout size that will approximately Nyquist
            # sample the map
            
            fwhm = {'030':0.54,
                    '044':0.45,
                    '070':0.22,
                    '100':0.16,
                    '143':0.12,
                    '217':0.08,
                    '353':0.08,
                    '545':0.08,
                    '857':0.08}[planck_freq]
                    
            pixel_size = fwhm/3
            n_pix = int(2*radius.value/pixel_size)
            
            pixels = str(n_pix)+','+str(n_pix)
            
            if verbose:
                print('Beginning Planck '+planck_freq)
                
            if not os.path.exists(galaxy + "/Planck/"):
                os.mkdir(galaxy + "/Planck/")

            if not os.path.exists(galaxy + "/Planck/" + planck_freq):
                os.mkdir(galaxy + "/Planck/" + planck_freq)
            
            if not os.path.exists(galaxy + "/Planck/" + planck_freq+'/outputs'):
                os.mkdir(galaxy + "/Planck/" + planck_freq+'/outputs')
                
            if 1 in steps:

                if verbose:
                    print("Downloading data")
    
                # Use SkyView to query and download a cutout around the
                # galaxy.
                
                planck_hdu = SkyView.get_images(position=galaxy,
                                                survey='Planck '+planck_freq+' I',
                                                radius=2*radius,
                                                pixels=pixels,
                                                cache=False,
                                                )
                
                planck_hdu[0].writeto(galaxy + "/Planck/" + planck_freq+'/outputs/'+galaxy+'.fits',
                                      overwrite=True)
                    
            if 2 in steps:
                 
                if verbose:
                    print('Converting to Jy')
             
                # Convert to Jy.
                 
                convert_to_jy(galaxy + "/Planck/" + planck_freq + "/outputs/"+galaxy+".fits",
                              planck_freq,
                              hdu_out=galaxy + "/Planck/"+galaxy+"_" + planck_freq+".fits")
                
        if original_radius is None:
            radius = None
        else:
            radius = original_radius.copy()
             
def convert_to_jy(hdu_in,frequency,hdu_out=None):
     
    """Convert from Planck K to Jy/pixel.
     
    Planck maps are units related to the CMB temperature. These can be
    converted to MJy/sr using conversion factors either in
    https://wiki.cosmos.esa.int/planckpla2015/index.php/UC_CC_Tables#LFI_Unit_Conversion_Tables
    for LFI, or Planck Collaboration IX (2013) for the HFI.
     
    Args:
        hdu_in (str or astropy.io.fits.PrimaryHDU): File name of Planck 
            .fits file, or an Astropy PrimaryHDU instance (i.e. the result 
            of ``fits.open(file)[0]``).
        frequency (str): Either '030', '044', '070', '100', '143', '217', 
            '353', '545', or '857'.
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
    
    # Convert from K to MJy/sr.
    
    k_conversion = {'030':23.5099,
                    '044':55.7349,
                    '070':129.1869,
                    '100':244.1,
                    '143':371.74,
                    '217':483.690,
                    '353':287.450,
                    '545':58.04,
                    '857':2.27}[frequency]
                    
    data *= k_conversion
     
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