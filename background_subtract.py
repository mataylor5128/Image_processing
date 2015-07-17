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

def weight_im(weight_image,image):
	for chip_num in range(im_nchip):
		image[0][chip_num] = image[0][chip_num]*weight_image[0][chip_num]/np.max(weight_image[0][chip_num])
	return image

def find_nearest(array,nvals,val):
	sorted = np.sort(np.abs(array))
	keep_vals = sorted[0:nvals]
	inds = []
	for value in keep_vals:
		inds.append(list(np.abs(array)).index(value))
	inds = np.array(inds)
	return inds

def median_subtract(back_cube,med_data,image,chip_num,q):
	for ii in range(im_size[0]):
		for jj in range(im_size[1]):
			back_val = np.median([back_cube[kk][chip_num][ii][jj] for kk in range(len(back_cube))])
			if np.isnan(back_val) == False:
				back_val_safe = back_val
				med_data[ii][jj] = back_val
			else:
				med_data[ii][jj] = back_val_safe

	image = image-med_data
	q.put([image, med_data])

def median_global_subtract(back_cube,med_data,image,chip_num,q):
	med_data[:][:] = np.median([back_cube[kk][chip_num] for kk in range(len(back_cube))])
	image = image-med_data
	q.put([image, med_data])

import pyfits
import numpy as np
import glob
from multiprocessing import Process, Queue
import matplotlib.pyplot as plt
# from scipy import interpolate
# from scipy import ndimage
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
import os,sys,subprocess
import time

start_time = time.time()

if len(sys.argv) < 2 :
	print "ERROR - missing subtraction method"
	exit()
if sys.argv[1] != "median":
	if sys.argv[1] != "median global":
		print "ERROR - background subtraction method must be one of [median, median global]"
		exit()

# im_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
# im_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
# do_tile = '3'
# do_dither = '3'
# do_filt = 'i'
# n_backs = 7

im_dir = sys.argv[2]
do_tile = sys.argv[3]
do_dither = sys.argv[4]
do_filt = sys.argv[5]
n_backs = int(sys.argv[6])
ss_im_out = sys.argv[7]
sky_im_out = sys.argv[8]

#back_dir = '/Volumes/MyPassport/masks/'
#mask_dir = '/Volumes/MyPassport/masks/'
#back_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
#mask_dir = '/Volumes/Q6/matt/2014A-0610/background_test_files/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
#back_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
#mask_dir = '/Volumes/Q6/matt/2014A-0610/pipeline/images/' #'/Volumes/Q6/matt/2014A-0610/pipeline/images/'
back_dir = im_dir
mask_dir = im_dir

test_image = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.fits'
test_weight = im_dir+'survey_t'+do_tile+'_d'+do_dither+'_'+do_filt+'_short.WEIGHT.fits'

im_h = pyfits.open(test_image)[0].header

#Open image file, and scale it by the weight map. ???SHOULD THE IMAGE BE SCALED????
print "\nProcessing image %s for background subtraction..." % test_image
print "\nLoading image..."
hdulist = pyfits.open(test_image)

im_nchip=0
for hdu in hdulist:
	if hdu.header['NAXIS'] != 0:
		hdu = np.array(hdu)
		im_dim= (4094, 2046) #hdu.shape
		im_ndim=len(im_dim)
		if im_ndim==2:
			if im_nchip==0: im_size=im_dim
			im_nchip+=1
hdulist.close()

im_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
#w_im_data = np.zeros((1,im_nchip,im_size[0],im_size[1]))
im_data = load_im(test_image,im_data)
#w_im_data =load_im(test_weight,w_im_data)

#Weight scale the original image
# im_data = weight_im(w_im_data,im_data)

#del w_im_data

# print np.max(im_data[0][0])

####################################################################
####Open background and mask files and scale them by weight maps####
####################################################################

#Determine the background files closest in time to the image
# temp_back_ims = glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.fits')
# temp_back_weights = glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.WEIGHT.fits')
# temp_mask_ims = glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.MASK.fits')
if do_dither == '1':
	temp_back_ims = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.fits')))
	temp_back_weights = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.WEIGHT.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.WEIGHT.fits')))
	temp_mask_ims = np.concatenate((glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.MASK.fits'),glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.MASK.fits')))
if do_dither != '1' and int(do_dither) <= 5:
	temp_back_ims = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.fits'),glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.fits')))
	temp_back_weights = np.concatenate((glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.WEIGHT.fits'),glob.glob(back_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.WEIGHT.fits'),glob.glob(back_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.WEIGHT.fits')))
	temp_mask_ims = np.concatenate((glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)-1)+'_*'+do_filt+'_*short.MASK.fits'),glob.glob(mask_dir+'survey_t*_*d'+do_dither+'_*'+do_filt+'_*short.MASK.fits'),glob.glob(mask_dir+'survey_t*_*d'+str(int(do_dither)+1)+'_*'+do_filt+'_*short.MASK.fits')))
if int(do_dither) > 5:
	temp_back_ims = glob.glob(back_dir+'survey_t*_*d*_*'+do_filt+'_*short.fits')
	temp_back_weights = glob.glob(back_dir+'survey_t*_*d*_*'+do_filt+'_*short.WEIGHT.fits')
	temp_mask_ims = glob.glob(mask_dir+'survey_t*_*d*_*'+do_filt+'_*short.MASK.fits')

# print temp_back_ims
# print temp_back_weights
# print temp_mask_ims

im_mjd = im_h['MJD-OBS']
print "\nImage observation date = %f" % im_mjd

#Delete tiles with bright, extended objects in them
#SCABS: tile1 (CenA), tile11 (wCen), tile24(wCen)
if do_dither == '1':
	gv_dithers = [do_dither,str(int(do_dither)+1)]
if do_dither != '1' and int(do_dither) <= 5:
	gv_dithers = [str(int(do_dither)-1),do_dither,str(int(do_dither)+1)]
if int(do_dither) > 5:
	gv_dithers = [str(ii+1) for ii in range(15)]
	
#print gv_dithers
bv_tiles = ['1','11','24']
for gv_d in gv_dithers:
	for bv in bv_tiles:
		if back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.fits' in temp_back_ims:
			temp_back_ims = np.delete(temp_back_ims, list(temp_back_ims).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.fits'),0)
			temp_back_weights = np.delete(temp_back_weights, list(temp_back_weights).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.WEIGHT.fits'), 0)
			temp_mask_ims = np.delete(temp_mask_ims, list(temp_mask_ims).index(back_dir+'survey_t'+bv+'_d'+gv_d+'_'+do_filt+'_short.MASK.fits'), 0)

for imfile in temp_back_ims:
	if imfile.replace('short.fits','short.MASK.fits') not in temp_mask_ims:
		temp_back_ims = np.delete(temp_back_ims, list(temp_back_ims).index(imfile), 0)		
for wfile in temp_back_weights:
	if wfile.replace('short.WEIGHT.fits','short.MASK.fits') not in temp_mask_ims:
		temp_back_weights = np.delete(temp_back_weights, list(temp_back_weights).index(wfile), 0)		
		
# print len(temp_back_ims), len(temp_back_weights), len(temp_mask_ims)
# print temp_back_ims
# print temp_back_weights
# print temp_mask_ims
#exit()
# print temp_back_ims
# print temp_back_weights

back_mjds = []
for ii in range(len(temp_back_ims)):
	back_mjds.append(float(subprocess.Popen("dfits "+temp_back_ims[ii]+" | fitsort MJD-OBS | awk 'NR>1 {print $2}'", shell=True, stdout=subprocess.PIPE).stdout.read().replace('\n','')))
back_mjds = np.array(back_mjds)

time_diffs = back_mjds - im_mjd
#print len(temp_back_ims), len(temp_back_weights), len(temp_mask_ims), len(time_diffs)
#print time_diffs

back_ims = []
back_weights = []
mask_ims = []

# print len(temp_back_ims), len(temp_back_weights), len(temp_mask_ims)

for idx in find_nearest(time_diffs,n_backs,0):
	back_ims.append(temp_back_ims[idx])
	back_weights.append(temp_back_weights[idx])
	mask_ims.append(temp_mask_ims[idx])

back_ims = np.array(back_ims)
back_weights = np.array(back_weights)
mask_ims = np.array(mask_ims)

# print back_ims
# print back_weights
# print mask_ims

print "\nUsing the following frames for the background subtraction: " 
print back_ims

back_mjds = []
for ii in range(n_backs):
	back_mjds.append(float(subprocess.Popen("dfits "+temp_back_ims[ii]+" | fitsort MJD-OBS | awk 'NR>1 {print $2}'", shell=True, stdout=subprocess.PIPE).stdout.read().replace('\n','')))

#Load background images and weight maps
back_im_cube = np.zeros( (len(back_ims), im_nchip, im_size[0], im_size[1]) )
for ii in range(len(back_ims)):
	print "\nLoading background image %i/%i into data cube..." % (ii+1,len(back_ims))
	back_im_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
	back_im_cube[ii] = load_im(back_ims[ii],back_im_temp)

#print back_im_cube
#print len(back_im_cube), len(back_im_cube[0]), len(back_im_cube[0][0]), len(back_im_cube[0][0][0])

#back_w_cube = np.zeros( (len(back_weights), im_nchip, im_size[0], im_size[1]) )

#Weight scale the backgroundi mages
# for ii in range(len(back_ims)):
# 	print "\nWeighting background image %i/%i..." % (ii+1,len(back_weights))
# 	back_weight_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
# 	back_w_data = load_im(back_weights[ii],back_weight_temp)
# 	back_im_cube[ii] = weight_im(back_w_data[0],back_im_cube[ii])
# 
# del back_w_data

#Load mask images and mask the background images
for ii in range(len(mask_ims)):
	print "\nMasking background image %i/%i..." % (ii+1,len(mask_ims))
	mask_im_temp = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )
	mask_im_data = load_im(mask_ims[ii],mask_im_temp)

	for jj in range(len(back_im_cube[0])):
		bv_mask = (mask_im_data[0][jj] == 1.)
		back_im_cube[ii][jj][bv_mask] = np.nan
	
del mask_im_data

########################################################################
####Calculate the background level from all of the background images####
########################################################################

# for ii in range(10):
#  	print "Chip #%i before subtraction:" % (ii+1)
# # 	print back_im_cube[0][ii][0][0], back_im_cube[1][ii][0][0], back_im_cube[2][ii][0][0]
# # 	print np.median([back_im_cube[0][ii][0][0], back_im_cube[1][ii][0][0], back_im_cube[2][ii][0][0]])
#  	print np.median([back_im_cube[kk][ii][:][:] for kk in range(len(back_im_cube))])
#  	print im_data[0][ii][0][0]

q = Queue()
results = np.zeros( (im_nchip, 2, im_size[0], im_size[1]) )
med_back_data = np.zeros( (1, im_nchip, im_size[0], im_size[1]) )

n_procs_max = 10
chip_start = 0
max_chips = 60
procs = []

while chip_start < max_chips:#im_nchips:
	print "Processing Chip #%i-%i..." % (chip_start+1, chip_start+1+n_procs_max)
	for chip_num in range(n_procs_max):
#  	  print "Processing Chip #%i" % (chip_num+1+chip_start)
	  if sys.argv[1] == "median":
		  procs.append(Process(target=median_subtract, args=(back_im_cube,med_back_data[0][chip_num+chip_start],im_data[0][chip_num+chip_start],chip_num+chip_start,q)))
	  if sys.argv[1] == "median global":
  		  procs.append(Process(target=median_global_subtract, args=(back_im_cube,med_back_data[0][chip_num+chip_start],im_data[0][chip_num+chip_start],chip_num+chip_start,q)))
	  #print "Proc %i start" % (chip_num+1+chip_start)
	  procs[chip_num+chip_start].start()
	  time.sleep(3.0)
	for chip_num in range(n_procs_max):
	  #print "Queue %i start" % (chip_num+1+chip_start)
	  results[chip_num+chip_start] = q.get()#result1 = q.get()
	  im_data[0][chip_num+chip_start] = results[chip_num+chip_start][0] #result1[0]
	  med_back_data[0][chip_num+chip_start] = results[chip_num+chip_start][1] #result1[1]
	chip_start+=n_procs_max

#Save the background subtracted image and the background image
hdulist = pyfits.open(test_image)

print "Saving image data..."
hdulist_out = hdulist
for ii in range(im_nchip):
	hdulist_out[ii+1].data = im_data[0][ii]
#hdulist_out.writeto('/Users/matt/Desktop/deleteme_image.fits',clobber=True)
hdulist_out.writeto('%s' % ss_im_out,clobber=True)
hdulist_out.close()

print "Saving background data..."
hdulist_back_out = hdulist
for ii in range(im_nchip):
	hdulist_back_out[ii+1].data = med_back_data[0][ii]
#hdulist_back_out.writeto('/Users/matt/Desktop/deleteme_back.fits',clobber=True)
hdulist_back_out.writeto('%s' % sky_im_out,clobber=True)
hdulist_back_out.close()

print "Done in %5.2f seconds." % (time.time() - start_time)

#Spline surface background subtraction
#Fit a spline to the surface of each chip, for each background image

