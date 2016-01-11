def load_im(image,array):
	jj=0
	hdulist = pyfits.open(image)
	for hdu in hdulist:
	  if hdu.header['NAXIS'] != 0:
	    im_dim=(4094, 2096) #hdu.shape
	    im_ndim=len(im_dim)
	    if im_ndim==2:
		  array[0,jj,:,:]=hdu.data
		  jj+=1
	hdulist.close()
	return array

import pyfits
from astropy.convolution import convolve, Gaussian2DKernel
import sys
import numpy as np

# im_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
im_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
back_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
mask_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'

do_tile = '1'
do_dither = '1'
do_filt = 'i'
n_backs = 3

back_dir = im_dir
mask_dir = im_dir

test_image = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.fits'
test_weight = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.WEIGHT.fits'
#test_mask = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.MASK.fits'
test_mask = im_dir+'sky_survey_t1_d1_g_short.004.fits'

im_h = pyfits.open(test_image)[0].header
chip_names = []

#Open image file, and scale it by the weight map. ???SHOULD THE IMAGE BE SCALED????
print "\nProcessing image %s for background subtraction..." % test_image
print "\nLoading image..."
hdulist = pyfits.open(test_image)

im_nchip=0
for hdu in hdulist:
	if hdu.header['NAXIS'] != 0:
#		hdu = np.array(hdu)
		im_dim = np.array(hdu.data).shape#(4096,2048)#(4094, 2046) #hdu.shape
		im_ndim=len(im_dim)
		chip_names.append(hdu.header['DETPOS'])
		if im_ndim==2:
			if im_nchip==0: im_size=im_dim
			im_nchip+=1
#hdulist.close()

# print chip_names.index('N21')
# exit()

im_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
im_data = load_im(test_image,im_data)
im_mask_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
im_masked = np.zeros((1,im_nchip,im_size[0],im_size[1]))

# if do_tile != '1':
# 	im_masked[:] = im_data[:]
im_mask_data = load_im(test_mask,im_mask_data)

gauss_kern = Gaussian2DKernel(2)
print gauss_kern.array
print "Convolving..."
im_mask_data[0][50] = convolve(im_mask_data[0][50], gauss_kern)

hdulist_back_out = hdulist
for ii in range(im_nchip):
      hdulist_back_out[ii+1].data = im_mask_data[0][ii]
hdulist_back_out.writeto('/Users/matt/Desktop/deleteme_mask_conv.fits',clobber=True)
hdulist_back_out.close()
exit()


# print "\nCreating masked image..."
# for ii in range(len(im_data[0])):
#  	bv_mask = (im_mask_data[0][ii] == 1.)
# 	im_masked[0][ii][bv_mask] = np.nan
