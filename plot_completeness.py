def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
    
def cutme(tab,cut,val):
	if cut == 'gtr':
		mask = np.where(tab > val,1, 0)
	elif cut == 'geq':
		mask = np.where(tab >= val,1, 0)
	elif cut == 'less':
		mask = np.where(tab < val,1, 0)
	elif cut == 'leq':
		mask = np.where(tab <= val,1, 0)
	elif cut == 'eq':
		mask = np.where(tab == val,1, 0)
	else:
		print "ERROR: no valid cutting option."
		exit()
	return mask
# 	tab = np.compress(mask,tab)
# 	return tab
	

def gen_fake_mags(mags_in,n):
	probs = np.logspace(1,1e-5,len(mags_in))
	probs = probs/np.sum(probs)
	return np.random.choice(np.sort(mags_in),size=n,replace=True,p=probs)

def mag_select(tab,mag_faint,mag_bright):
	tab = np.compress(cutme(tab,'less',mag_faint),tab)
	tab = np.compress(cutme(tab,'geq',mag_bright),tab)
	return tab

def comp_calc(mag_chk,dmag):
	N_real = []
	N_rec = []
	mags = []
	while mag_chk < np.max(mag_real):
		mags.append(np.mean([mag_chk,mag_chk+dmag]))
# 		print np.mean([mag_chk,mag_chk+dmag]), len(mag_select(tbdata1.field('MAG_PSF'),mag_chk+dmag,mag_chk)), len(mag_select(fake_mags,mag_chk+dmag,mag_chk))
		N_real.append(float(len(mag_select(tbdata1.field('MAG_PSF'),mag_chk+dmag,mag_chk))))
		N_rec.append(float(len(mag_select(fake_mags,mag_chk+dmag,mag_chk))))
		mag_chk += dmag
	N_real = np.array(N_real)
	N_rec = np.array(N_rec)
	probs = N_rec/N_real
	probs = probs/np.max(probs)#(N_rec-N_real)/N_real*100.
	return np.array(mags), probs

from astropy.io import fits
from astropy.coordinates import ICRS, match_coordinates_sky#SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import random

#Load the mock catalogue with all fake sources
file1 = '/Volumes/MyPassport/pipeline_v4/survey_tile1_i_psf.ldac'
#Load output catalogue with list of fake sources that were actually detected
file2 = '/Volumes/MyPassport/pipeline_v4/survey_tile1_u_psf_ALIGNi.ldac'

hdu1 = fits.open(file1)
hdu2 = fits.open(file2)

tbdata1 = hdu1[2].data
tbdata2 = hdu2[2].data


mag_real = np.compress(cutme(tbdata1.field('MAG_PSF'),'less',99.),tbdata1.field('MAG_PSF'))
ra_real  = np.compress(cutme(tbdata1.field('MAG_PSF'),'less',99.),tbdata1.field('ALPHA_J2000'))
dec_real = np.compress(cutme(tbdata1.field('MAG_PSF'),'less',99.),tbdata1.field('DELTA_J2000'))

n = len(mag_real)#200000
fake_mags = gen_fake_mags(mag_real,n)

mags, probs = comp_calc(np.min(mag_real),0.1)

print "50 percent completeness = %.2f" % (mags[find_nearest(probs,0.5)])

plt.figure()
plt.plot([np.min(mags),np.max(mags)],[0,0],'k--')
plt.plot(mags,probs,'o')
plt.xlim(np.min(mags),np.max(mags))
plt.ylim(0,1.1)

dcoord = 500 #arcsec
dcoord = dcoord/3600.
# ra_chk = np.max(ra_real)
dec_chk = np.min(dec_real)
x = np.linspace(np.max(ra_real),np.min(ra_real),int((np.max(ra_real)-np.min(ra_real))/dcoord)-1)
y = np.linspace(np.min(dec_real),np.max(dec_real),int((np.max(dec_real)-np.min(dec_real))/dcoord)-1)
xplot, yplot = np.meshgrid(x,y,sparse=True)
coord_dats = np.zeros((int((np.max(ra_real)-np.min(ra_real))/dcoord)+1, int((np.max(dec_real)-np.min(dec_real))/dcoord+1)))
# dec_dats = np.zeroes( (np.max(ra_real)-np.min(ra_real)/dcoord), (np.max(ra_real)-np.min(ra_real)/dcoord) )
# xx = 0
yy = 0
dmag=0.05
zz = []
while dec_chk < np.max(dec_real):
	xx = 0
	ra_chk = np.max(ra_real)
	print dec_chk, np.max(dec_real)
	mask = cutme(dec_real,'geq',dec_chk)
	mag_temp = np.compress(mask,mag_real)
	ra_temp  = np.compress(mask,ra_real)
	dec_temp  = np.compress(mask,dec_real)
	print "Check 1 " + str(len(mag_temp))
	mask = cutme(dec_temp,'less',dec_chk+dcoord)
	mag_temp = np.compress(mask,mag_temp)
	ra_temp  = np.compress(mask,ra_temp)
	dec_temp  = np.compress(mask,dec_temp)
	print "Check 2 " + str(len(mag_temp))
	while ra_chk > np.min(ra_real):
		print ra_chk, np.min(ra_real)
		mask = cutme(ra_temp,'leq',ra_chk)
		mag_temp = np.compress(mask,mag_temp)
		ra_temp  = np.compress(mask,ra_temp)
		dec_temp  = np.compress(mask,dec_temp)
		print "Check 3 " + str(len(mag_temp))
		mask = cutme(ra_temp,'gtr',ra_chk-dcoord)
		mag_temp = np.compress(mask,mag_temp)
		ra_temp  = np.compress(mask,ra_temp)
		dec_temp  = np.compress(mask,dec_temp)
		print "Check 4 " + str(len(mag_temp))

		n = len(mag_temp)
		print n
		if n > 0:
			fake_temp = gen_fake_mags(mag_temp,n)
			mags, probs = comp_calc(np.min(mag_temp),dmag)
			print len(mags), "50 percent comp limit at (%.4f,%.4f) = %.2f" % (np.mean([ra_chk,ra_chk+dcoord]),np.mean([dec_chk,dec_chk+dcoord]),mags[find_nearest(probs,0.5)])
			coord_dats[xx][yy] = mags[find_nearest(probs,0.5)]
			
		elif n == 0:
			coord_dats[xx][yy] = None#np.nan
		xx += 1
		ra_chk -= dcoord
	yy += 1
	dec_chk += dcoord

plt.figure()
#plt.imshow(coord_dats)
# plt.contourf(x,y,coord_dats)
plt.colorbar()
plt.show()
exit()
# ra1 = tbdata1.field('ALPHA_J2000')
# dec1 = tbdata1.field('DELTA_J2000')

c1 = ICRS(ra=tbdata1.field('ALPHA_J2000'), dec=tbdata1.field('DELTA_J2000'), unit=(u.degree, u.degree))
c2 = ICRS(ra=tbdata2.field('ALPHA_J2000'), dec=tbdata2.field('DELTA_J2000'), unit=(u.degree, u.degree))

idx, d2d, d3d = match_coordinates_sky(c1,c2)

matches = c2[idx]
dra = (matches.ra - c1.ra).arcsec
ddec = (matches.dec - c1.dec).arcsec
rad = np.sqrt(dra**2+ddec**2)

mask = np.where(rad >= 1.0,1,0)
mat_ra = np.compress(mask,tbdata1.field('ALPHA_J2000'))
mat_dec = np.compress(mask,tbdata1.field('DELTA_J2000'))



# plt.figure()
# plt.plot(mat1.field('ALPHA_J2000'),mat1.field('DELTA_J2000'),'bo',ms=3)
# plt.plot(mat2.field('ALPHA_J2000'),mat2.field('DELTA_J2000'),'ro',ms=3)
# plt.show()