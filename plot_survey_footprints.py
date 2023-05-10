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
cmass_cat_file1 = '/Users/jackelvinpoole/SDSS/DR12/BOSS/LSS/galaxy_DR12v5_CMASS_North.fits.gz'
cmass_cat_file2 = '/Users/jackelvinpoole/SDSS/DR12/BOSS/LSS/galaxy_DR12v5_CMASS_South.fits.gz'
desi_tile_file = '/Users/jackelvinpoole/DESI/desi_tiles.fits'
lsst_footprint_file = '/Users/jackelvinpoole/DESC/footprint/opsim_nvisits_g.fits.gz'

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

lsst_hpix_ring = fio.read(lsst_footprint_file)['I'].flatten()
lsst_nside = hp.npix2nside(len(lsst_hpix_ring))
lsst_hpix_nest = lsst_hpix_ring[hp.nest2ring(lsst_nside, np.arange(len(lsst_hpix_ring )))]
lsst = hsp.HealSparseMap(nside_coverage=nside_coverage, healpix_map=lsst_hpix_nest )
lsst = lsst.degrade(nside_sparse,reduction='max')
lsst_ra,lsst_dec = lsst.valid_pixels_pos(lonlat=True)
lsst_ra[lsst_ra > 180.] -= 360.
lsst_ra *= -1.
lsst_ra = lsst_ra*np.pi/180.
lsst_dec = lsst_dec*np.pi/180.


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


plt.figure('phot',figsize=(10,5))
plt.subplot(111, projection="mollweide")
lsst_color = 'wheat'
plt.scatter(lsst_ra,lsst_dec,marker='.',color=lsst_color,s=4.0,label='LSST')
plt.scatter(des_ra,des_dec,marker='.',color='r',s=4.0,label='DES')
plt.scatter(unions_ra,unions_dec,marker='.',color='c',s=4.0,label='UNIONS')
plt.scatter(hsc_ra,hsc_dec,marker='.',color='orange',s=1.,label='HSC')
plt.scatter(kids_ra,kids_dec,marker='.',color='b',s=4.0,label='KiDS-1000')
lsstbox = dict(boxstyle='round', facecolor='white', edgecolor=lsst_color, lw=3.)
redbox = dict(boxstyle='round', facecolor='white', edgecolor='r', lw=3.)
cyanbox = dict(boxstyle='round', facecolor='white', edgecolor='c', lw=3.)
orangebox = dict(boxstyle='round', facecolor='white', edgecolor='orange', lw=3.)
bluebox = dict(boxstyle='round', facecolor='white', edgecolor='b', lw=3.)
plt.text(130.*np.pi/180.,-45.*np.pi/180.,'Rubin-LSST',fontsize=16, bbox=lsstbox)
plt.text(-140.*np.pi/180.,-45.*np.pi/180.,'DES',fontsize=16, bbox=redbox)
plt.text(-30.*np.pi/180.,60.*np.pi/180.,'UNIONS',fontsize=16, bbox=cyanbox)
plt.text(120.*np.pi/180.,20.*np.pi/180.,'HSC',fontsize=16, bbox=orangebox)
plt.text(50.*np.pi/180.,-25.*np.pi/180.,'KiDS-1000',fontsize=16, bbox=bluebox)
plt.savefig('phot_future.png')
plt.close()

#SPEC-z
cmass_fits1 = fio.FITS(cmass_cat_file1)
cmass_ra1  = cmass_fits1[-1]['RA'].read()
cmass_dec1 = cmass_fits1[-1]['Dec'].read()
cmass_fits2 = fio.FITS(cmass_cat_file2)
cmass_ra2  = cmass_fits2[-1]['RA'].read()
cmass_dec2 = cmass_fits2[-1]['Dec'].read()
cmass_ra = np.append(cmass_ra1,cmass_ra2)
cmass_dec = np.append(cmass_dec1,cmass_dec2)
cmass_ra[cmass_ra > 180.] -= 360.
cmass_ra *= -1.
cmass_ra  = cmass_ra*np.pi/180.
cmass_dec = cmass_dec*np.pi/180.


desi_fits = fio.FITS(desi_tile_file)
desi_ra  = desi_fits[-1]['RA'].read()
desi_dec = desi_fits[-1]['Dec'].read()
desi_ra[desi_ra > 180.] -= 360.
desi_ra *= -1.
desi_ra  = desi_ra*np.pi/180.
desi_dec = desi_dec*np.pi/180.
desi_ra = desi_ra + np.random.rand(len(desi_ra))*0.1 - 0.05
desi_dec = desi_dec + np.random.rand(len(desi_dec))*0.1 - 0.05


plt.figure('spec',figsize=(10,5))
plt.subplot(111, projection="mollweide")
plt.scatter(desi_ra,desi_dec,marker='s',color='k',s=70.0,label='DESI')
plt.scatter(cmass_ra,cmass_dec,marker='s',color='g',s=40.0,label='BOSS',alpha=0.3)

"""
desi1 = np.loadtxt('desi_polygon_1',unpack=True,delimiter=',')
desi_line1 = matplotlib.lines.Line2D(desi1[0], desi1[1], lw=3., color='k')
desi2 = np.loadtxt('desi_polygon_2',unpack=True,delimiter=',')
desi_line2 = matplotlib.lines.Line2D(desi2[0], desi2[1], lw=3., color='k')
desi3 = np.loadtxt('desi_polygon_3',unpack=True,delimiter=',')
desi_line3 = matplotlib.lines.Line2D(desi3[0], desi3[1], lw=3., color='k')

ax = plt.gca()
ax.add_line(desi_line1)
ax.add_line(desi_line2)
ax.add_line(desi_line3)
"""

greenbox = dict(boxstyle='round', facecolor='white', edgecolor='g', lw=3.)
blackbox = dict(boxstyle='round', facecolor='white', edgecolor='k',lw=3.)

plt.text(-190.*np.pi/180.,30.*np.pi/180.,'BOSS',fontsize=16, bbox=greenbox)
plt.text(70.*np.pi/180.,-30.*np.pi/180.,'DESI',fontsize=16, bbox=blackbox)

plt.savefig('specz_current.png')
plt.close()



