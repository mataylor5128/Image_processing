# Create RGB images for big astronomical images (DECam,VISTA,VST)
# Copyright 2014, 2015 Roberto Pablo Munoz, All rights reserved
#
# This product includes software developed by
# Roberto Pablo Mu~noz in Santiago, Chile
# (https://github.com/rpmunoz/DECam).

import sys,os
import pyfits
import astropy
import aplpy
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import skimage.io as io
#from skimage import exposure
#from PIL import Image
#from PIL import ImageEnhance

stack_dir='/Volumes/Q6/matt/stacks'
#stack_dir='/Users/matt/Desktop'
weight_dir='/Volumes/Q6/matt/stacks'
stack_tile=['2']#,'2','3']
stack_tile_ref='2'
stack_filter=['i','r','g']
stack_filter_ref='i'
stack_version='1'

#stack_im_file=['survey_tileX_%s_short_ALIGN%s.fits' % (stack_filter[0],stack_filter_ref),'survey_tileX_%s_short.fits' % stack_filter[1],'survey_tileX_%s_short_ALIGN%s.fits' % (stack_filter[2],stack_filter_ref)]
#stack_im_file=['ss_survey_tileX_%s_short.004.fits' % (stack_filter[0]),'ss_survey_tileX_%s_short_ALIGN%s.004.fits' % (stack_filter[1],stack_filter_ref),'ss_survey_tileX_%s_short_ALIGN%s.004.fits' % (stack_filter[2],stack_filter_ref)]
stack_im_file=['survey_tileX_%s_short.fits' % (stack_filter[0]),'survey_tileX_%s_short_ALIGN%s.fits' % (stack_filter[1],stack_filter_ref),'survey_tileX_%s_short_ALIGN%s.fits' % (stack_filter[2],stack_filter_ref)]
#stack_weight_file=['survey_tileX_%s_short_ALIGN%s.WEIGHT.fits' % (stack_filter[0],stack_filter_ref),'survey_tileX_%s_short.WEIGHT.fits' % stack_filter[1],'survey_tileX_%s_short_ALIGN%s.WEIGHT.fits' % (stack_filter[2],stack_filter_ref)]
#stack_weight_file=['ss_survey_tileX_%s_short.004.WEIGHT.fits' % (stack_filter[0]),'ss_survey_tileX_%s_short_ALIGN%s.004.WEIGHT.fits' % (stack_filter[1],stack_filter_ref),'ss_survey_tileX_%s_short_ALIGN%s.004.WEIGHT.fits' % (stack_filter[2],stack_filter_ref)]
stack_weight_file=['survey_tileX_%s_short.WEIGHT.fits' % (stack_filter[0]),'survey_tileX_%s_short_ALIGN%s.WEIGHT.fits' % (stack_filter[1],stack_filter_ref),'survey_tileX_%s_short_ALIGN%s.WEIGHT.fits' % (stack_filter[2],stack_filter_ref)]
#stack_mask_file=['survey_tileX_%s_short_ALIGN%s.MASK.fits' % (stack_filter[0],stack_filter_ref),'survey_tileX_%s_short.MASK.fits' % stack_filter[1],'survey_tileX_%s_short_ALIGN%s.MASK.fits' % (stack_filter[2],stack_filter_ref)]
#stack_mask_file=['ss_survey_tileX_%s_short.004.MASK.fits' % stack_filter[0],'ss_survey_tileX_%s_short.004.MASK.fits' % stack_filter[1],'survey_tileX_%s_short.004.MASK.fits' % stack_filter[2]]
stack_mask_file=['survey_tileX_%s_short.MASK.fits' % stack_filter[0],'survey_tileX_%s_short.MASK.fits' % stack_filter[1],'survey_tileX_%s_short.MASK.fits' % stack_filter[2]]
# stack_im_file=['survey_tileX_%s_short.fits' % stack_filter[0],'survey_tileX_%s_short_ALIGN%s.fits' % (stack_filter[1],stack_filter_ref),'survey_tileX_%s_short_ALIGN%s.fits' % (stack_filter[2],stack_filter_ref)]
# stack_weight_file=['survey_tileX_%s_short.WEIGHT.fits' % stack_filter[0],'survey_tileX_%s_short_ALIGN%s.WEIGHT.fits' % (stack_filter[1],stack_filter_ref),'survey_tileX_%s_short_ALIGN%s.WEIGHT.fits' % (stack_filter[2],stack_filter_ref)]
# #stack_mask_file=['survey_tileX_%s_short_ALIGN%s.MASK.fits' % (stack_filter[0],stack_filter_ref),'survey_tileX_%s_short.MASK.fits' % stack_filter[1],'survey_tileX_%s_short_ALIGN%s.MASK.fits' % (stack_filter[2],stack_filter_ref)]
# stack_mask_file=['survey_tileX_%s_short.MASK.fits' % stack_filter[0],'survey_tileX_%s_short.MASK.fits' % stack_filter[1],'survey_tileX_%s_short.MASK.fits' % stack_filter[2]]

# stack_im_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_im_file] for im_tile in stack_tile]
# stack_weight_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_weight_file] for im_tile in stack_tile]
# stack_mask_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_mask_file] for im_tile in stack_tile]
# stack_rgb_file_full=[stack_dir+'/stack_tile'+im_tile+'_rgb.fits' for im_tile in stack_tile]
stack_im_file_full = [[stack_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_im_file] for im_tile in stack_tile]
stack_weight_file_full = [[weight_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_weight_file] for im_tile in stack_tile]
stack_mask_file_full = [[weight_dir+'/'+im_file.replace('tileX','tile'+im_tile) for im_file in stack_mask_file] for im_tile in stack_tile]
stack_rgb_file_full=[stack_dir+'/stack_tile'+im_tile+'_rgb.fits' for im_tile in stack_tile]

print stack_im_file_full
print stack_weight_file_full
print stack_mask_file_full
print stack_rgb_file_full
<<<<<<< Updated upstream

# for ii in range(len(stack_im_file_full[0])):
# 	if os.path.exists(stack_im_file_full[0][ii]) != True:
#  		print stack_im_file_full[0][ii]+" does not exist."
#  		exit()
# 	if os.path.exists(stack_weight_file_full[0][ii]) != True:
#  		print stack_im_file_full[0][ii]+" does not exist."
#  		exit()
# 	if os.path.exists(stack_mask_file_full[0][ii]) != True:
#  		print stack_im_file_full[0][ii]+" does not exist."
#  		exit()
=======
exit()
>>>>>>> Stashed changes

for i in range(len(stack_im_file_full)):

	stack_im_file=stack_im_file_full[i]
	stack_weight_file=stack_weight_file_full[i]
	stack_mask_file=stack_mask_file_full[i]
	stack_rgb_file=stack_rgb_file_full[i]
	stack_rgb_limit=np.zeros((3,2), dtype=np.float32)

	hist_nbin=200
	hist_percentile=[0.25,99.5] #[0.25,99.8] #[0.25,99.5]  #([0.25,99.5],[0.25,99.55],[0.22,99.8])
	
	hdulist = pyfits.open(stack_im_file[0])
	im_h=hdulist[0].header
	hdulist.close()
	
	nx = int(im_h['NAXIS1'])
	ny = int(im_h['NAXIS2'])
	im_data_cube = np.zeros((3, ny, nx), dtype=np.float32)
	
	if not os.path.exists(stack_rgb_file):
	
		print '\nCreating rgb cube ', stack_rgb_file
		for j in range(len(stack_filter)):
	
			print '\nProcessing image ', stack_im_file[j], ' - filter ', stack_filter[j]
			im_file=stack_im_file[j]
			weight_file=stack_weight_file[j]
			mask_file=stack_mask_file[j]
		
			if os.path.exists(im_file):
				print '\nReading image file ', im_file
				hdulist = pyfits.open(im_file)
				im_data=hdulist[0].data
				im_h=hdulist[0].header
				hdulist.close()
		
				print 'Reading weight file ', weight_file
				hdulist = pyfits.open(weight_file)
				weight_data=hdulist[0].data
				hdulist.close()

				print "Check1"
		
				bv_weight = (weight_data==0.)

				print "CHeck2"
				im_data[bv_weight]=np.nan
		
				print 'Image size along X-axis ', im_h['NAXIS1']
				print 'Image size along Y-axis ', im_h['NAXIS2']
		
				ypad = ny - im_data.shape[0]
				xpad = nx - im_data.shape[1]
				pad_width = ((0, ypad), (0, xpad))

				im_data_pad = np.pad(im_data, pad_width, mode='constant', constant_values=[np.nan])
				im_data_cube[j,:,:]=im_data_pad
			else:
				print 'The image file does not exists ', im_file
	
				nim_data = np.empty((ny, nx), dtype=np.float32)
				nim_data[:] = np.nan
				try:
					gv = stack_tile.index(stack_tile_ref)
					im_file=stack_im_file_full[gv][j]
					weight_file=stack_weight_file_full[gv][j]
					mask_file=stack_mask_file_full[gv][j]

					print '\nReading image file ', im_file
					hdulist = pyfits.open(im_file)
					im_data=hdulist[0].data
					hdulist.close()
		
					print 'Reading weight file ', weight_file
					hdulist = pyfits.open(weight_file)
					weight_data=hdulist[0].data
					hdulist.close()
		
					print 'Reading mask file ', mask_file
					hdulist = pyfits.open(mask_file)
					mask_data=hdulist[0].data
					hdulist.close()

					gv = np.logical_and(weight_data>0,mask_data==0)
			
					if np.count_nonzero(gv)>1:
						sky_data=im_data[gv]
						hist_limit= np.abs(np.percentile(sky_data, 0.001))
						hist_xmin=-1.*hist_limit
						hist_xmax=1.*hist_limit
						hist, bin_edges = np.histogram( sky_data, bins=np.linspace(hist_xmin, hist_xmax, hist_nbin), range=(hist_xmin, hist_xmax))
						plt.bar(bin_edges[:-1], hist, width = bin_edges[1]-bin_edges[0], linewidth=0., log=True)
						plt.xlim(hist_xmin, hist_xmax)
						plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
						plt.xlabel('Flux')
						plt.ylabel('Number count')
						plt.savefig(stack_dir+'/stack_tile'+stack_tile_ref+'_histogram_sky.pdf')

						weight_file=stack_weight_file[0]
						print 'Reading weight file ', weight_file
						hdulist = pyfits.open(weight_file)
						weight_data=hdulist[0].data
						hdulist.close()

						print 'Random sample from the reference image'
						gv=(weight_data>0)
						nim_data[gv]=np.random.choice(sky_data, np.count_nonzero(gv))
						im_data_cube[j,:,:]=nim_data

				except ValueError:
					print 'Adding gaussian model'
					im_data_cube[j,:,:]= 6. * np.random.randn(ny, nx)
					im_data_cube[j,:,:][bv_weight]=np.nan

		
		print 'Saving the rgb cube file ', stack_rgb_file
		pyfits.writeto(stack_rgb_file, im_data_cube, header=im_h)
	

	im_rgb_file=stack_rgb_file.replace('.fits','_asinh_v'+stack_version+'.jpg')
	if not os.path.exists(im_rgb_file):

		print "Reading the RGB fits cube file: ", stack_rgb_file		
		hdulist = pyfits.open(stack_rgb_file)
		im_data_cube=hdulist[0].data
		im_h_cube=hdulist[0].header
		hdulist.close()
	
		w, h = 2*matplotlib.figure.figaspect(1.2)
		fig = plt.figure(figsize=(w,h))
		for j in range(len(stack_filter)):
			print '\nProcessing filter ', stack_filter[j]
			im_data=im_data_cube[j,:,:]
			gv=np.isfinite(im_data)

			if np.count_nonzero(gv)>1:

				if os.path.exists(stack_im_file[j]):
					print 'Computing the rgb limit using the actual data'

					stack_rgb_limit[j,:]= np.percentile(im_data[gv], hist_percentile)  #[plot_xmin, plot_xmax]
					if stack_rgb_limit[j,0]==stack_rgb_limit[j,1]:
						print 'The data seems to contain a constant value. Filter '+stack_filter[j]
						continue

					hist_xmin=np.amin(im_data[gv])
					hist_xmax=np.amax(im_data[gv])
	
					print 'Percentiles for computing RGB min and max ', hist_percentile
					print 'Histogram min and max ', hist_xmin, hist_xmax
					print 'RGB image min and max ', stack_rgb_limit[j,:]
	
					hist, bin_edges = np.histogram( im_data, bins=np.linspace(hist_xmin, hist_xmax, hist_nbin), range=(hist_xmin, hist_xmax))
		
					plt.subplot(3,2,1+2*j)
					plt.bar(bin_edges[:-1], hist, width = bin_edges[1]-bin_edges[0], linewidth=0., log=True)
					plt.xlim(hist_xmin, hist_xmax)
					plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
					plt.title('Tile '+stack_tile[i]+' - '+stack_filter[j])
					plt.xlabel('Flux')
					plt.ylabel('Number count')
	
					if stack_filter[j]=='i': #always the first filter in list
						hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], hist_nbin), range=stack_rgb_limit[j,:] )
						stack_rgb_ratio_ref= [ hist[0]*1./np.amax(hist), hist[-1]*1./np.amax(hist) ]
						gv_max=np.argmax(hist)
						stack_rgb_limit_ratio_ref= (bin_edges[gv_max]-stack_rgb_limit[j,0])*1./(stack_rgb_limit[j,1]-stack_rgb_limit[j,0])
						print 'stack_rgb_ratio_ref ', stack_rgb_ratio_ref
						print 'stack_rgb_limit_ratio_ref ', stack_rgb_limit_ratio_ref
					else:
						stack_rgb_ratio=stack_rgb_ratio_ref
						print 'stack_rgb_ratio ', stack_rgb_ratio
	
						delta=stack_rgb_limit[j,1]-stack_rgb_limit[j,0]
						if stack_filter[j]=='u':
							stack_rgb_limit[j,0] -= delta/8.
							stack_rgb_limit[j,1] += 2*delta
						else:
							stack_rgb_limit[j,0] -= delta/8.
							stack_rgb_limit[j,1] += delta/4.
	
						print "Guess RGB image min and max ", stack_rgb_limit[j,:]
				
						hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], 2*hist_nbin), range=stack_rgb_limit[j,:] )
						gv_max=np.argmax(hist)
						gv_limit_max=gv_max + np.argmin( np.abs(hist[gv_max:-1]*1./np.amax(hist) - stack_rgb_ratio[1]) )
						gv_limit_min=np.argmin( np.abs( (bin_edges[gv_max]- bin_edges[0:gv_max])*1./(bin_edges[gv_limit_max]-bin_edges[0:gv_max]) - stack_rgb_limit_ratio_ref) )
						stack_rgb_limit[j,:]=bin_edges[[gv_limit_min,gv_limit_max]]
						print "Computed RGB limit ", stack_rgb_limit[j,:]

				else:
					print 'Computing the rgb limit using the reference tile'

					gv_tile = stack_tile.index(stack_tile_ref)
					gv_filter = stack_filter.index(stack_filter_ref)
					hdulist = pyfits.open(stack_rgb_file_full[gv_tile])
					im_h=hdulist[0].header
					hdulist.close()

					stack_rgb_limit[j,:]=np.asarray(im_h['RGB_'+stack_filter[j]].split(',')).astype(np.float)
					delta=stack_rgb_limit[j,1]-stack_rgb_limit[j,0]
					stack_rgb_limit[j,0] -= delta/4.
					stack_rgb_limit[j,1] += delta/4.

					hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], hist_nbin), range=stack_rgb_limit[j,:] )
					gv_max=np.argmax(hist)

					stack_rgb_limit[j,1]= np.asarray(im_h['RGB_'+stack_filter[j]].split(',')).astype(np.float)[1]/np.asarray(im_h['RGB_'+stack_filter_ref].split(',')).astype(np.float)[1]*stack_rgb_limit[gv_filter,1]
					stack_rgb_limit[j,0]= bin_edges[ np.argmin( np.abs( (bin_edges[gv_max]- bin_edges[0:gv_max])*1./(stack_rgb_limit[j,1]-bin_edges[0:gv_max]) - stack_rgb_limit_ratio_ref) ) ]
					print "Computed RGB limit ", stack_rgb_limit[j,:]
					
				hist, bin_edges = np.histogram( im_data, bins=np.linspace(stack_rgb_limit[j,0], stack_rgb_limit[j,1], hist_nbin), range=stack_rgb_limit[j,:] )
	
				print 'Adding RGB limits to header'
				print stack_rgb_limit[j,:]
				im_h_cube['RGB_'+stack_filter[j]]= ','.join(np.char.mod('%.2f', stack_rgb_limit[j,:]))

				plt.subplot(3,2,1+2*j+1)
				plt.bar(bin_edges[:-1], hist, width = bin_edges[1]-bin_edges[0], linewidth=0., log=True)
				plt.xlim(stack_rgb_limit[j,0], stack_rgb_limit[j,1])
				plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
				plt.title('Tile '+stack_tile[i]+' - '+stack_filter[j])
				plt.xlabel('Flux')
				plt.ylabel('Number count')

			else:
				print '\nNo valid pixels for filter '+stack_filter[j]

		print 'Updating the rgb cube file ', stack_rgb_file
		pyfits.writeto(stack_rgb_file, im_data_cube, header=im_h_cube, clobber=True)
	
		plt.tight_layout()
		plt.savefig(stack_dir+'/stack_tile'+stack_tile[i]+'_histogram_v'+stack_version+'.pdf')
		plt.close(fig)
	
		print "stack_rgb_limit: \n", stack_rgb_limit
		
		im_rgb_file=stack_rgb_file.replace('.fits','_asinh_v'+stack_version+'.jpg')
		print 'Creating RGB image file ', im_rgb_file
		aplpy.make_rgb_image(stack_rgb_file, im_rgb_file, stretch_r='arcsinh', stretch_g='arcsinh', stretch_b='arcsinh', vmin_r=stack_rgb_limit[0,0], vmin_g=stack_rgb_limit[1,0], vmin_b=stack_rgb_limit[2,0], vmax_r=stack_rgb_limit[0,1], vmax_g=stack_rgb_limit[1,1], vmax_b=stack_rgb_limit[2,1], vmid_r=-0.07, vmid_g=-0.07, vmid_b=-0.07, make_nans_transparent=True, embed_avm_tags=True)

		im_rgb_file=stack_rgb_file.replace('.fits','_linear_v'+stack_version+'.jpg')
		print 'Creating RGB image file ', im_rgb_file
		aplpy.make_rgb_image(stack_rgb_file, im_rgb_file, stretch_r='linear', stretch_g='linear', stretch_b='linear', vmin_r=stack_rgb_limit[0,0], vmin_g=stack_rgb_limit[1,0], vmin_b=stack_rgb_limit[2,0], vmax_r=stack_rgb_limit[0,1], vmax_g=stack_rgb_limit[1,1], vmax_b=stack_rgb_limit[2,1], make_nans_transparent=True, embed_avm_tags=True)

