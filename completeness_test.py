def add_source(n,x,y,prev_src,xwidth):	
	dmax = 2*np.sqrt(2*max_d**2)
	print "Adding sources..."
	sources = []
	ii = 0
	while ii < n:
		#Select random coordinate in chip's view
		x_src = int(random.uniform(x,x+xwidth))
		y_src = int(random.uniform(y,y+2249))
		#Check segmentation map to see if anything is nearby
		if np.max(seg_image[0].section[y_src-max_d:y_src+max_d,x_src-max_d:x_src+max_d]) == 0:
			#If this isn't the first chip, check that the new source doesn't
			#come close to source on previous chip (e.g. sources near edges.
			if len(prev_src) > 0:
				for key in prev_src:
					D = []
					for kk in range(len(source_coords[key])):
						xpix = source_coords[key][kk][0]
						ypix = source_coords[key][kk][1]
						D.append(np.sqrt((xpix-x_src)**2+(ypix-y_src)**2))
					if np.min(D) > dmax:
						#If it's not the first source on-chip, then check if new
						#source is space-invading previous source.
						if len(sources) > 0:
							D = []
							for jj in range(len(sources)):
								D.append(np.sqrt((x_src-sources[jj][0])**2+(y_src-sources[jj][1])**2))
							if np.min(D) > dmax:
									sources.append([x_src,y_src])
									ii = ii + 1
							#Okay, nothing close to previous chip sources. Is this the first source on chip?
						if len(sources) == 0:
							sources.append([x_src,y_src])
							ii = ii + 1

			#Is this the first chip?
			if len(prev_src) == 0:
				if len(sources) > 0:
					D = []
 					for jj in range(len(sources)):
 						D.append(np.sqrt((x_src-sources[jj][0])**2+(y_src-sources[jj][1])**2))
					if np.min(D) > dmax:
						sources.append([x_src,y_src])
						ii = ii + 1
				#Is this the first artificial source?
				if len(sources) == 0:
					sources.append([x_src,y_src])
					ii = ii + 1
	return sources
	
def cut_stamp(x,y,type):
	if type == 'image':
		stamp = image[0].section[y:y+2249,x:x+4314]
		seg_stamp = seg_image[0].section[y:y+2249,x:x+4314]
	if type == 'source':
		stamp = source_image[0].data[y:y+2249,x:x+4314]
		seg_stamp = seg_image[0].section[y:y+2249,x:x+4314]
	return stamp, seg_stamp

def gauss_source(x,y):
	point_x = np.arange(x-max_d,x+max_d)
	point_y = np.arange(y-max_d,y+max_d)
	point_x,point_y = np.meshgrid(point_x,point_y)

	#amp = random.uniform(5e2,2e5)
	amp, mag = get_cnts()
	gauss = Gaussian2D(amp, x, y, 0.5, 0.5)

	data_2D = gauss(point_x,point_y) 
	#data_2D = data_2D + 1000 * (np.random.rand(2*max_d, 2*max_d) - 0.5)
	gauss_kernel = Gaussian2DKernel(3)

	smoothed_data_gauss = convolve(data_2D, gauss_kernel)
	#smoothed_data_gauss = smoothed_data_gauss + 100 * (np.random.rand(2*max_d, 2*max_d) - 0.5)

	return smoothed_data_gauss, mag

def moffat_source(x,y):
	point_x = np.arange(x-max_d,x+max_d)
	point_y = np.arange(y-max_d,y+max_d)
	point_x,point_y = np.meshgrid(point_x,point_y)

	#amp = random.uniform(5e2,2e5)
	amp, mag = get_cnts()
	
	moffat = Moffat2D(amp, x, y, 0.5, 0.5)

	data_2D = moffat(point_x,point_y) 
	#data_2D = data_2D + 1000 * (np.random.rand(2*max_d, 2*max_d) - 0.5)
	moffat_kernel = Model2DKernel(Moffat2D,3,3)

	smoothed_data_moffat = convolve(data_2D, moffat_kernel)

	return smoothed_data_moffat, mag

def get_cnts():
	t_exp = 100.
	gain = 4.5
	jy = 3730
	dlamb = 0.14
	mag = random.uniform(15.,25.)
	cnts = (10**(-0.4*mag)*jy)*(1.51e7)*dlamb*np.pi*2.**2 #phot/s
	cnts = cnts*t_exp/gain #adu
	
	return cnts, mag

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import random
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.modeling.models import Gaussian2D
#from astropy.modeling.functional_models import Moffat2D
from PIL import Image

#Load Image
image = pyfits.open('/Volumes/MyPassport/pipeline_v2/full_runs/2MASS/ngc5128_tile1_g.fits')
source_image = image
header = image[0].header
seg_image = pyfits.open('/Volumes/MyPassport/pipeline_v2/full_runs/2MASS/ngc5128_tile1_g.CHECK_SEGMENTATION.fits')
seg_header = seg_image[0].header

source_coords = {}
n = 2500
max_d = 12

#Row 1: 2 chips
source_coords['S31'] = add_source(n,9270,25622,source_coords,4314)
source_coords['S29'] = add_source(n,17775,25645,source_coords,4314)

#Row 2: 4 chips
x1 = 7150
x2 = 24221
dx = int((x2-x1)/4)
source_coords['S28'] = add_source(n,x1,23367,source_coords,dx)
source_coords['S27'] = add_source(n,x1+dx,23379,source_coords,dx)
source_coords['S26'] = add_source(n,x1+2*dx,23391,source_coords,dx)
source_coords['S25'] = add_source(n,x1+3*dx,23404,source_coords,dx)

#Row 3: 5 chips
x1 = 5030
x2 = 26350
dx = int((x2-x1)/5))
source_coords['S24'] = add_source(n,x1,21112,source_coords,dx)
source_coords['S23'] = add_source(n,x1+dx,21123,source_coords,dx)
source_coords['S22'] = add_source(n,x1+2*dx,21134,source_coords,dx)
source_coords['S21'] = add_source(n,x1+3*dx,21145,source_coords,dx)
source_coords['S20'] = add_source(n,x1+4*dx,21156,source_coords,dx)

#Row 4: 6 chips
x1 = 2909
x2 = 28480
dx = int((x2-x1)/6)
source_coords['S19'] = add_source(n,x1,18754,source_coords,dx)
source_coords['S18'] = add_source(n,x1+dx,18770,source_coords,dx)
source_coords['S17'] = add_source(n,x1+2*dx,18786,source_coords,dx)
source_coords['S16'] = add_source(n,x1+3*dx,18802,source_coords,dx)
source_coords['S15'] = add_source(n,x1+4*dx,18818,source_coords,dx)
source_coords['S14'] = add_source(n,x1+5*dx,18834,source_coords,dx)

#Row 5: 6 chips
x1 = 2914
x2 = 28483
dx = int((x2-x1)/6)
source_coords['S13'] = add_source(n,x1,16607,source_coords,dx)
source_coords['S12'] = add_source(n,x1+dx,16617,source_coords,dx)
source_coords['S11'] = add_source(n,x1+2*dx,16628,source_coords,dx)
source_coords['S10'] = add_source(n,x1+3*dx,16638,source_coords,dx)
source_coords['S09'] = add_source(n,x1+4*dx,16649,source_coords,dx)
source_coords['S08'] = add_source(n,x1+5*dx,16660,source_coords,dx)

#Row 6: 7 chips
x1 = 822
x2 = 30661
dx = int((x1-x2)/2)
source_coords['S07'] = add_source(n,x1,14134,source_coords,dx)
source_coords['S06'] = add_source(n,x1+dx,14170,source_coords,dx)
source_coords['S05'] = add_source(n,x1+2*dx,14207,source_coords,dx)
source_coords['S04'] = add_source(n,x1+3*dx,14233,source_coords,dx)
source_coords['S03'] = add_source(n,x1+4*dx,14134,source_coords,dx)
source_coords['S02'] = add_source(n,x1+5*dx,14134,source_coords,dx)
source_coords['S01'] = add_source(n,x1+6*dx,14354,source_coords,dx)

fileout = open('source_coords_mags.txt','w')
print >> fileout, '#xpix		ypix	mag'

#mags = []
for key in source_coords:
	print key
	for ii in range(len(source_coords[key])):
		xpix = source_coords[key][ii][0]
		ypix = source_coords[key][ii][1]
		bkgrnd = source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d]
#		source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d] = moffat_source(xpix,ypix) + bkgrnd
		source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d], mag = gauss_source(xpix,ypix)
		source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d] = source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d] + bkgrnd
#		source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d] = source_image[0].data[ypix-max_d:ypix+max_d,xpix-max_d:xpix+max_d] + bkgrnd
		print >> fileout, '%i		%i		%6.3f' % (xpix, ypix, mag)

fileout.close()
pyfits.writeto('source_image.fits',source_image[0].data,header,clobber=True)
exit()

chips = {}
#Row 1: 2 chips 
S31, seg_S31 = cut_stamp(9270,25622,'image')
chips['S31.fits'] = S31
S29, seg_S29 = cut_stamp(17775,25645,'image')
chips['S29.fits'] = S29
S28, seg_S28 = cut_stamp(7150,23367,'image')
chips['S28.fits'] = S28
S27, seg_S27 = cut_stamp(11464,23379,'image')
chips['S27.fits'] = S27
S26, seg_S26 = cut_stamp(15778,23391,'image')
chips['S26.fits'] = S26
S25, seg_S25 = cut_stamp(19906,23404,'image')
chips['S25.fits'] = S25

source_chips = {}
src_S31, seg_src_31 = cut_stamp(9270,25622,'source')
source_chips['src_S31.fits'] = src_S31
src_S29, seg_src_29 = cut_stamp(17775,25645,'source')
source_chips['src_S29.fits'] = src_S29
src_S28, src_seg_S28 = cut_stamp(7150,23367,'source')
source_chips['S28.fits'] = S28
src_S27, src_seg_S27 = cut_stamp(11464,23379,'source')
source_chips['S27.fits'] = S27
src_S26, src_seg_S26 = cut_stamp(15778,23391,'source')
source_chips['S26.fits'] = S26
src_S25, src_seg_S25 = cut_stamp(19906,23404,'source')
source_chips['S25.fits'] = S25

for key, value in sorted(chips.items()):
	pyfits.writeto(key,value,header,clobber=True)

for key, value in sorted(source_chips.items()):
	pyfits.writeto(key,value,header,clobber=True)

exit()

##############
command='sex '+sex_stack_im_file+' -c sex_config/ctio_decam_stack.sex -CATALOG_NAME '+sex_stack_cat_file+' -WEIGHT_IMAGE '+sex_stack_weight_file+' -MAG_ZEROPOINT 30. -XML_NAME '+sex_stack_xml_file+' -CHECKIMAGE_TYPE '+sex_stack_checkimage_type+' -CHECKIMAGE_NAME '+sex_stack_checkimage_file
##############


chk = 0
while chk ==0:
	x, y, result = add_source(9270,25622)
	if result == "Source Nearby":
		chk = 1
	print x, y, result



exit()

#Slice image into regions
chips = {}
seg_chips = {}
#Row 1: 2 chips 
S31, seg_S31 = cut_stamp(9270,25622)
chips['S31.fits'] = S31
seg_chips['seg_S31.fits'] = seg_S31
S29, seg_S29 = cut_stamp(17775,25645)
chips['S29.fits'] = S29
seg_chips['seg_S29.fits'] = seg_S29
#Row 2: 4 chips
S28, seg_S28 = cut_stamp(7150,23367)
chips['S28.fits'] = S28
seg_chips['seg_S28.fits'] = seg_S28
S27, seg_S27 = cut_stamp(11464,23367)
chips['S27.fits'] = S27
seg_chips['seg_S27.fits'] = seg_S27
S26, seg_S26 = cut_stamp(15778,23367)
chips['S26.fits'] = S26
seg_chips['seg_S26.fits'] = seg_S26
S25, seg_S25 = cut_stamp(19906,23404)
chips['S25.fits'] = S25
seg_chips['seg_S25.fits'] = seg_S25

for key, value in sorted(chips.items()):
	pyfits.writeto(key,value,header,clobber=True)
for key, value in sorted(chips.items()):
	pyfits.writeto(key,value,header,clobber=True)

exit()
#Randomly add source locations, using segmentation map to avoid crowded regions

#Randomly assign magnitudes to the artificial sources

#Use psfex-modelled psfs to add accurately modelled artificial stellar sources to image

#Run SEXZTRACTOR on 


# Generate fake data
# gauss = Gaussian2D(1, 0, 0, 1, 1)
# 
# x = np.arange(-100,101)
# y = np.arange(-100,101)
# x, y = np.meshgrid(x,y)
# 
# data_2D = gauss(x, y) + 0.1 * (np.random.rand(201, 201) - 0.5)
# 
# gauss_kernel = Gaussian2DKernel(5)
# 
# smoothed_data_gauss = convolve(data_2D, gauss_kernel)
# 
# plt.figure()
# plt.imshow(data_2D, vmin=-0.01, vmax=0.08, origin='lower', interpolation='None')
# 
# plt.colorbar()
# 
# plt.figure()
# plt.imshow(smoothed_data_gauss, vmin=-0.01, vmax=0.08, origin='lower', interpolation='None')
# 
# plt.colorbar()
# 
# 
# plt.show()
# exit()


#Load image
image = pyfits.open('/Volumes/MyPassport/pipeline_v2/distortion_check/2MASS/full_run/images/ngc5128_tile2_i.fits')

#print image[0].header['CRVAL1'], image[0].header['CRVAL2'], image[0].header['NAXIS1']
print image[0].header['CRPIX1'], image[0].header['CRPIX2'], image[0].header['NAXIS1']
#print image[0].data[int(image[0].header['CRPIX1']),int(image[0].header['CRPIX2'])]


min_ra = image[0].header['CRVAL1'] - (np.abs(image[0].header['CD1_1'])*image[0].header['CRPIX1'])*(1./np.cos(image[0].header['CRVAL2']*np.pi/180.))
max_ra = image[0].header['CRVAL1'] + (np.abs(image[0].header['CD1_1'])*image[0].header['CRPIX1'])*(1./np.cos(image[0].header['CRVAL2']*np.pi/180.))
print min_ra, max_ra
min_dec = image[0].header['CRVAL2'] - (np.abs(image[0].header['CD2_2'])*image[0].header['CRPIX2'])
max_dec = image[0].header['CRVAL2'] + (np.abs(image[0].header['CD2_2'])*image[0].header['CRPIX2'])
print min_dec, max_dec

min_x = 1.0
max_x = 2.*image[0].header['CRPIX1']
min_y = 1.0
max_y = 2.*image[0].header['CRPIX2']
print min_x, max_x
print min_y, max_y

#Randomly select point
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

x = np.arange(-1000,1001)
y = np.arange(-1000,1001)
x,y = np.meshgrid(x,y)


point_coords = []
final_2D = []
final_smooth = []
for ii in range(1):

	point = [random.uniform(min_x,max_x),random.uniform(min_y,max_y)]
	print point
	point_coords.append(point)

	point_x = np.arange(np.min(point[0])-100,np.max(point[0])+101)
	point_y = np.arange(np.min(point[1])-100,np.max(point[1])+101)
	point_x,point_y = np.meshgrid(point_x,point_y)

	amp = random.uniform(1,100)
	gauss = Gaussian2D(amp, point[0], point[1], 0.5, 0.5)

	data_2D = gauss(point_x,point_y)
	gauss_kernel = Gaussian2DKernel(5)

	smoothed_data_gauss = convolve(data_2D, gauss_kernel)

	if ii == 0:
		final_2D = data_2D
		final_smooth = smoothed_data_gauss
		im1 = ax1.imshow(final_2D, vmin=-0.01, vmax=0.08, origin='lower')#, interpolation='None')

		im2 = ax2.imshow(final_smooth, vmin=-0.01, vmax=0.08, origin='lower')#, interpolation='None')

		fig1.colorbar(im1)
		fig2.colorbar(im2)
	
		plt.show()
		exit()


	if ii >= 0:
		final_2D = np.concatenate([final_2D,data_2D])
		final_smooth = np.concatenate([final_smooth,smoothed_data_gauss])
#	final_2D.append(data_2D)

#	ax1.imshow(data_2D, vmin=-0.01, vmax=0.08, origin='lower', interpolation='None')
im1 = ax1.imshow(final_2D, vmin=-0.01, vmax=0.08, origin='lower')#, interpolation='None')

im2 = ax2.imshow(final_smooth, vmin=-0.01, vmax=0.08, origin='lower')#, interpolation='None')

fig1.colorbar(im1)
fig2.colorbar(im2)

plt.show()
exit()
