import numpy as np
import matplotlib.pyplot as plt 
import fitsio as fio 
import healpy as hp 
import healsparse as hsp

nside_coverage = 32 
nside_sparse = 128

#stage 4 lensing
des_mask_file = '/Users/jackelvinpoole/DES/cats/y3/masks/MASK_Y3LSSBAOSOF_22_3_v2p2_y3_format.fits.gz'
kids_cat_file = '/Users/jackelvinpoole/KIDS/KiDS_DR4.zhQAYSRk.1_ugriZYJHKs_SOM_gold_WL_cat.fits'
kids_cat_file = '/Users/jackelvinpoole/KIDS/KiDS_radec_first100000.fits'
hsc_cat_file = '/Users/jackelvinpoole/HSC/random_radec_100000_wide.fits'
unions_cat_file = '/Users/jackelvinpoole/UNIONS/DataReleases/v1.0/ShapePipe/unions_shapepipe_2022_v1.0.fits'
boss_cat_file1 = '/Users/jackelvinpoole/SDSS/DR12/BOSS/LSS/galaxy_DR12v5_CMASS_North.fits.gz'
boss_cat_file2 = '/Users/jackelvinpoole/SDSS/DR12/BOSS/LSS/galaxy_DR12v5_CMASS_South.fits.gz'


#### load maps

des_hpix_4096 = fio.read(des_mask_file)['HPIX']
hp_des = np.zeros(hp.nside2npix(4096))+hp.UNSEEN
hp_des[hp.ring2nest(4096,des_hpix_4096)] = 1. #just set it to 1 for this, these are only approx
des = hsp.HealSparseMap(nside_coverage=nside_coverage, healpix_map=hp_des)
des = des.degrade(nside_sparse,reduction='max')
des_ra,des_dec = des.valid_pixels_pos(lonlat=True)
des_ra[des_ra > 180.] -= 360.
des_ra *= -1.
des_ra = des_ra*np.pi/180.
des_dec = des_dec*np.pi/180.


kids_fits = fio.FITS(kids_cat_file)
#kids_ra  = kids_fits[-1]['RAJ2000'].read()
#kids_dec = kids_fits[-1]['DECJ2000'].read()
#subsample here if using full cat
kids_ra  = kids_fits[-1]['s_ra'].read()
kids_dec = kids_fits[-1]['s_dec'].read()
kids_ra[kids_ra > 180.] = kids_ra[kids_ra > 180.] - 360.
kids_ra *= -1.
kids_ra  = kids_ra*np.pi/180.
kids_dec = kids_dec*np.pi/180.

hsc_fits = fio.FITS(hsc_cat_file)
hsc_ra  = hsc_fits[-1]['ra'].read()
hsc_dec = hsc_fits[-1]['dec'].read()
hsc_ra[hsc_ra > 180.] -= 360.
hsc_ra *= -1.
hsc_ra  = hsc_ra*np.pi/180.
hsc_dec = hsc_dec*np.pi/180.

unions_fits = fio.FITS(unions_cat_file)
unions_ra  = unions_fits[-1]['RA'].read()
unions_dec = unions_fits[-1]['Dec'].read()
nunions=len(unions_ra)
#select = np.random.choice(np.arange(nunions),size=100000, replace=False)
select = np.random.randint(0,nunions,size=100000)
unions_ra = unions_ra[select]
unions_dec = unions_dec[select]
unions_ra[unions_ra > 180.] -= 360.
unions_ra *= -1.
unions_ra  = unions_ra*np.pi/180.
unions_dec = unions_dec*np.pi/180.


plt.figure('phot',figsize=(10,5))
plt.subplot(111, projection="mollweide")
plt.scatter(des_ra,des_dec,marker='.',color='r',s=4.0,label='DES')
plt.scatter(unions_ra,unions_dec,marker='.',color='c',s=4.0,label='UNIONS')
plt.scatter(hsc_ra,hsc_dec,marker='.',color='orange',s=1.,label='HSC')
plt.scatter(kids_ra,kids_dec,marker='.',color='b',s=4.0,label='KiDS-1000')
redbox = dict(boxstyle='round', facecolor='white', edgecolor='r', lw=3.)
cyanbox = dict(boxstyle='round', facecolor='white', edgecolor='c', lw=3.)
orangebox = dict(boxstyle='round', facecolor='white', edgecolor='orange', lw=3.)
bluebox = dict(boxstyle='round', facecolor='white', edgecolor='b', lw=3.)
plt.text(-140.*np.pi/180.,-45.*np.pi/180.,'DES',fontsize=16, bbox=redbox)
plt.text(-30.*np.pi/180.,60.*np.pi/180.,'UNIONS',fontsize=16, bbox=cyanbox)
plt.text(120.*np.pi/180.,20.*np.pi/180.,'HSC',fontsize=16, bbox=orangebox)
plt.text(50.*np.pi/180.,-25.*np.pi/180.,'KiDS-1000',fontsize=16, bbox=bluebox)
plt.savefig('phot_current.png')
plt.close()




