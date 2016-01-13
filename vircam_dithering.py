def draw_chips(x,y,chip_x0,gap_lng0):	
	#Draw out the 16 VIRCam chips and return them as a collections of patches
	#(x0,y0) = RA, Dec of the bottom corner of each chip, such that (0,0) is the bottom-left corner of array
	#Input x,y are coordinates of the centre of the array
	chips = []
	for ii in range(4):
		for jj in range(4):
			y0 = (y+gap_sht/2.+(jj-2.)*chip_y+(jj-2.)*gap_sht)
			chip_x = chip_x0/np.cos(y0*np.pi/180.)
			gap_lng = gap_lng0/np.cos(y0*np.pi/180.)
			x0 = (x+gap_lng/2.+(ii-2.)*chip_x+(ii-2.)*gap_lng)
			chips.append(patches.Rectangle([x0,y0],chip_x,chip_y))
			hits_ra.append([x0,(x0+chip_x)])
			hits_dec.append([y0,y0+chip_y])
	return PatchCollection(chips,alpha=0.1)

def get_tiles():
	#Return the set of coordinates corresponding to the centres of each tile
	#Need to cover 10 tiles in sets of 2
	#Columns of 3 1/2 FoV to the East, and 1/2 FoV to the West
	#Column of 4 running 1.5 FoV N, 0.5 N, 0.5 S, 1.5 S
	#			t8
	#   t2				t5
	#			t1
	#	t7				t6
	#			t3
	#	t9				t4
	#			t10
	#Pair as: (t1, t2); (t3,t4); (t5,t6); (t7,t8); (t9,t10)
	#(t1,t2,t7);(t3,t6,t4);(t5,t8);(t9,t10)
	tiles = [[cena_ra,cena_dec+fov_dec/2.], 									 		#Tile 1 coordinates 
			[cena_ra+fov_ra/np.cos((cena_dec+fov_dec)*np.pi/180.),cena_dec+fov_dec], 	#Tile 2 coordinates 
			[cena_ra,cena_dec+fov_dec/2.],												#Tile 3 coordinates 
			[cena_ra-fov_ra/np.cos((cena_dec-fov_dec)*np.pi/180.),cena_dec-fov_dec],	#Tile 4 coordinates 
			[cena_ra-fov_ra/np.cos((cena_dec+fov_dec)*np.pi/180.),cena_dec+fov_dec],	#Tile 5 coordinates 
			[cena_ra-fov_ra/np.cos((cena_dec)*np.pi/180.),cena_dec],					#Tile 6 coordinates 
			[cena_ra+fov_ra/np.cos((cena_dec)*np.pi/180.),cena_dec],					#Tile 7 coordinates 
			[cena_ra,cena_dec+1.5*fov_dec],												#Tile 8 coordinates 
			[cena_ra+fov_ra/np.cos((cena_dec-fov_dec)*np.pi/180.),cena_dec-fov_dec],	#Tile 9 coordinates 
			[cena_ra,cena_dec-1.5*fov_dec]												#Tile 10 coordinates 
				]
	return tiles

def get_spiral(coords):
	#Draw the Fermat spiral and return the coordinates
	theta = 137.508*np.pi/180.
	a = d/(2.*np.sqrt(n-1))
	b = 1.0 #Try 4-5 arcmin
	pointings = []
	for ii in range(int(n)):
		ra0 = a*np.sqrt(ii)*np.cos(b*ii*theta)
		dec0 = a*np.sqrt(ii)*np.sin(b*ii*theta)
		pointings.append([coords[0]+ra0,coords[1]+dec0])
	return pointings

def get_paw(spir1, spir2):
	#For each point in the two input Fermat spirals, cycle between them running through the 6-pointing pawprint
	coords_out = []
	dra = 0.95*chip_x
	ddec = 0.475*chip_y
	for ii in range(len(spir1)):
		dcoords1 = [[-0.5*dra/np.cos((spir1[ii][1]-ddec)*np.pi/180.),-1.*ddec],
					[-0.5*dra/np.cos(spir1[ii][1]*np.pi/180.)       ,0.],
					[-0.5*dra/np.cos((spir1[ii][1]+ddec)*np.pi/180.),ddec],
					[0.5*dra/np.cos((spir1[ii][1]+ddec)*np.pi/180.) ,ddec],
					[0.5*dra/np.cos(spir1[ii][1]*np.pi/180.)        ,0.],
					[0.5*dra/np.cos((spir1[ii][1]-ddec)*np.pi/180.),-1.*ddec]
					]
		dcoords2 = [[-0.5*dra/np.cos((spir2[ii][1]-ddec)*np.pi/180.),-1.*ddec],
					[-0.5*dra/np.cos(spir2[ii][1]*np.pi/180.)       ,0.],
					[-0.5*dra/np.cos((spir2[ii][1]+ddec)*np.pi/180.),ddec],
					[0.5*dra/np.cos((spir2[ii][1]+ddec)*np.pi/180.) ,ddec],
					[0.5*dra/np.cos(spir2[ii][1]*np.pi/180.)        ,0.],
					[0.5*dra/np.cos((spir2[ii][1]-ddec)*np.pi/180.),-1.*ddec]
					]
		for jj in range(6):
			coords_out.append([spir1[ii][0]+dcoords1[jj][0],spir1[ii][1]+dcoords1[jj][1]])
			coords_out.append([spir2[ii][0]+dcoords2[jj][0],spir2[ii][1]+dcoords2[jj][1]])
	return coords_out

# def draw_footprint(point_ra,point_dec,chip_x,gap_sht):
# 	if sys.argv[1] == 'rect':
# 		points_x = [point_ra,point_ra+dra,point_ra+dra,point_ra]
# 		points_y = [point_dec,point_dec,point_dec+ddec,point_dec+ddec]
# 		p1 = PatchCollection(draw_chips(point_ra,point_dec,chip_x,gap_sht),alpha=0.1)
# 		p2 = PatchCollection(draw_chips(dra,point_dec,chip_x,gap_sht),alpha=0.1)
# 		p3 = PatchCollection(draw_chips(dra,ddec,chip_x,gap_sht),alpha=0.1)
# 		p4 = PatchCollection(draw_chips(point_ra,ddec,chip_x,gap_sht),alpha=0.1)
# 		collection = [p1,p2,p3,p4]
# 
# 	if sys.argv[1] == 'rect_centre':
# 		points_x = [point_ra,point_ra-1.*dra,point_ra+dra,point_ra+dra,point_ra-1.*dra]
# 		points_y = [point_dec,point_dec-1.*ddec,point_dec-1.*ddec,point_dec+ddec,point_dec+ddec]
# 		p1 = PatchCollection(draw_chips(point_ra,point_dec,chip_x,gap_sht),alpha=0.05)
# 		p2 = PatchCollection(draw_chips(point_ra-dra,point_dec-ddec,chip_x,gap_sht),alpha=0.1)
# 		p3 = PatchCollection(draw_chips(point_ra+dra,point_dec-ddec,chip_x,gap_sht),alpha=0.1)
# 		p4 = PatchCollection(draw_chips(point_ra+dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
# 		p5 = PatchCollection(draw_chips(point_ra-dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
# 		collection = [p1,p2,p3,p4,p5]
# 
# 	if sys.argv[1] == 'hexagon':
# 		points_x = [point_ra,point_ra-1.*dra,point_ra-1.*dra,point_ra,dra,dra]
# 		points_y = [point_dec,point_dec+ddec,point_dec+2.*ddec,point_dec+3.*ddec,point_dec+2.*ddec,point_dec+ddec]
# 		p1 = PatchCollection(draw_chips(point_ra,point_dec,chip_x,gap_sht),alpha=0.1)
# 		p2 = PatchCollection(draw_chips(point_ra-dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
# 		p3 = PatchCollection(draw_chips(point_ra-dra,point_dec+2.*ddec,chip_x,gap_sht),alpha=0.1)
# 		p4 = PatchCollection(draw_chips(point_ra,point_dec+3.*ddec,chip_x,gap_sht),alpha=0.1)
# 		p5 = PatchCollection(draw_chips(point_ra+dra,point_dec+2.*ddec,chip_x,gap_sht),alpha=0.1)
# 		p6 = PatchCollection(draw_chips(point_ra+dra,point_dec+ddec,chip_x,gap_sht),alpha=0.1)
# 		collection = [p1,p2,p3,p4,p5,p6]
# 	
# 	return points_x, points_y, collection

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import sys

dit = 20.

cena_ra = (13.+25./60.+27.61507/3600.)*15. 
cena_dec = (-43.-1./60.-8.8053/3600.)

xlim = [cena_ra+7.5/np.cos(cena_dec*np.pi/180.),cena_ra-7.5/np.cos(cena_dec*np.pi/180.)]
ylim = [cena_dec-6.5,cena_dec+5.0]

wcen_ra = (13.+26./60.+47.28/3600.)*15. 
wcen_dec = -47-28./60.-46.1/3600.

gal_coords = [[(13.+27./60.+37.3/3600.)*15.,-41.-28./60.-50./3600.],
			  [(13.+21./60.+47.3/3600.)*15.,-45.-3./60.-42.0/3600.],
			  [(13.+37./60.+39.0/3600.)*15.,-42.-50./60.-49.0/3600.],
			  [(13.+13./60.+11.9/3600.)*15.,-43.-15./60.-56.0/3600.],
			  [(13.+13./60.+9.1/3600.)*15.,-44.-53./60.-24.0/3600.],
			  [(13.+34./60.+47.3/3600.)*15.,-45.-32./60.-51.0/3600.],
			  [(13.+45./60.+0.0/3600.)*15.,-41.-51./60.-40.0/3600.],
			  [(13.+37./60.+25.3/3600.)*15.,-39.-53./60.-48.0/3600.],
			  [(13.+10./60.+32.9/3600.)*15.,-46.-59./60.-27.0/3600.],
			  [(13.+5./60.+2.1/3600.)*15.,-40.-4./60.-58.0/3600.],
			  [(13.+30./60+14.26/3600.)*15.,-41.-53./60.-35.8/3600.], #Crnojevic
			  [(13.+29./60.+57.34/3600.)*15.,-41.-52./60.-22.6/3600.],#] #Crnojevic
			  [(13.+22./60.+02.10/3600)*15.,-42.-32./60.-05.9/3600.], #karachentsev
			  [(13.+22./60.+12.30/3600)*15.,-42.-44./60.-04.6/3600.], #karachentsev
			  [(13.+30./60+21.5/3600.)*15.,-42.-11./60.-33.0/3600.], #Crnojevic
			  [(13.+23./60+02.6/3600.)*15.,-41.-47./60.-10.0/3600.], #Crnojevic
			  [(13.+19./60+52.4/3600.)*15.,-41.-59./60.-37.0/3600.], #Crnojevic
			  [(13.+25./60+57.6/3600.)*15.,-41.-05./60.-39.0/3600.], #Crnojevic
			  [(13.+26./60+28.7/3600.)*15.,-43.-33./60.-24.0/3600.], #Crnojevic
			  [(13.+33./60+34.1/3600.)*15.,-41.-36./60.-28.0/3600.], #Crnojevic
			  [(13.+33./60+01.5/3600.)*15.,-42.-31./60.-48.0/3600.], #Crnojevic
			  ]
			  
crno_coords = [[(13.+24./60.+10.90/3600.)*15.,-42.-08./60.-24.5/3600.],
			 	[(13.+27./60.+11./3600)*15.,  -42.-29./60.-57./3600.],
			 	[(13.+26./60.+49.16/3600)*15.,-43.-00./60.-01./3600.],
			 	[(13.+27./60.+01.38/3600)*15.,-43.-03./60.-28./3600.],
			 	[(13.+24./60.+10.20/3600)*15.,-43.-23./60.-52./3600.],
			 	[(13.+22./60.+52.7/3600)*15., -43.-32./60.-45./3600.],
			 	[(13.+23./60.+45.4/3600)*15., -43.-39./60.-21.5/3600.],
#			 	[(13.+26./60.+28./3600)*15.,  -43.-33./60.-15.6/3600.],
			 	[(13.+27./60.+57.27/3600)*15.,-43.-33./60.-07.8/3600.],
			 	[(13.+28./60.+16.9/3600)*15., -43.-42./60.-18.4/3600.],
			 	[(13.+26./60.+07.42/3600)*15.,-43.-42./60.-28.9/3600.],
#			 	[(13.+22./60.+02.10/3600)*15.,-42.-32./60.-05.9/3600.],
# 			 	[(13.+22./60.+12.30/3600)*15.,-42.-44./60.-04.6/3600.],
			 	[(13.+23./60.+27.62/3600)*15.,-41.-49./60.-08.8/3600.],
			 	[(13.+24./60.+53.83/3600)*15.,-40.-45./60.-39.0/3600.],
			 	[(13.+23./60.+54.60/3600)*15.,-40.-49./60.-58.1/3600.],
			 	[(13.+22./60.+45.31/3600)*15.,-41.-32./60.-22.8/3600.],
			 	[(13.+21./60.+18.21/3600)*15.,-40.-52./60.-41.5/3600.],
			 	[(13.+25./60.+13.76/3600)*15.,-41.-29./60.-40.7/3600.],
			 	[(13.+28./60.+18.15/3600)*15.,-41.-26./60.-41.2/3600.],
#			 	[(13.+23./60.+02.41/3600)*15.,-41.-47./60.-05.3/3600.],
			 	[(13.+33./60.+01.45/3600)*15.,-41.-09./60.-04.0/3600.],
			 	[(13.+35./60.+23.10/3600)*15.,-41.-20./60.-41.7/3600.],
			 	[(13.+35./60.+30.76/3600)*15.,-41.-20./60.-34.5/3600.],
			 	[(13.+34./60.+04.94/3600)*15.,-41.-15./60.-16.4/3600.],
			 	[(13.+31./60.+33.93/3600)*15.,-41.-37./60.-47.7/3600.],
			 	[(13.+30./60.+04.28/3600)*15.,-41.-37./60.-46.2/3600.],
			 	[(13.+29./60.+45.10/3600)*15.,-41.-38./60.-49.2/3600.],
			 	[(13.+34./60.+13.94/3600)*15.,-41.-47./60.-03.9/3600.],
			 	[(13.+35./60.+03.05/3600)*15.,-41.-55./60.-06.8/3600.],
			 	[(13.+28./60.+34.78/3600)*15.,-41.-54./60.-13.2/3600.],
			 	[(13.+27./60.+44.87/3600)*15.,-42.-06./60.-13.0/3600.],
			 	[(13.+28./60.+04.47/3600)*15.,-42.-08./60.-33.9/3600.],
			 	[(13.+28./60.+50.47/3600)*15.,-42.-05./60.-04.5/3600.],
			 	[(13.+33./60.+09.75/3600)*15.,-42.-14./60.-19.6/3600.],
			 	[(13.+31./60.+45.60/3600)*15.,-42.-30./60.-38.2/3600.],
			 	[(13.+35./60.+38.04/3600)*15.,-42.-23./60.-57.2/3600.],
			 	[(13.+31./60.+20.11/3600)*15.,-42.-32./60.-12.2/3600.],
			 	[(13.+15./60.+24.80/3600)*15.,-41.-30./60.-16.6/3600.],
			 	[(13.+13./60.+36.54/3600)*15.,-41.-36./60.-10.2/3600.],
			 	[(13.+14./60.+44.81/3600)*15.,-41.-42./60.-28.5/3600.],
			 	[(13.+19./60.+36.36/3600)*15.,-41.-33./60.-01.5/3600.],
			 	[(13.+14./60.+44.71/3600)*15.,-41.-42./60.-29.1/3600.],
			 	[(13.+12./60.+45.12/3600)*15.,-41.-49./60.-57.1/3600.],
			 	[(13.+14./60.+02.70/3600)*15.,-42.-00./60.-07.1/3600.],
			 	[(13.+13./60.+58.90/3600)*15.,-41.-58./60.-03.4/3600.],
			 	[(13.+17./60.+14.47/3600)*15.,-41.-56./60.-32.1/3600.],
			 	[(13.+19./60.+21.47/3600)*15.,-42.-03./60.-41.1/3600.],
			 	[(13.+14./60.+08.08/3600)*15.,-42.-04./60.-08.9/3600.],
			 	[(13.+12./60.+22.28/3600)*15.,-42.-18./60.-43.4/3600.],
			 	[(13.+13./60.+34.16/3600)*15.,-42.-11./60.-09.9/3600.],
			 	[(13.+13./60.+36.48/3600)*15.,-42.-14./60.-07.3/3600.],
			 	[(13.+13./60.+59.85/3600)*15.,-42.-21./60.-38.9/3600.],
			 	[(13.+16./60.+42.35/3600)*15.,-42.-24./60.-05.4/3600.],
			 	[(13.+15./60.+03.05/3600)*15.,-42.-32./60.-18.0/3600.],
			 	[(13.+14./60.+21.93/3600)*15.,-42.-30./60.-40.4/3600.],
			 	[(13.+20./60.+23.32/3600)*15.,-42.-32./60.-08.9/3600.],
			 	[(13.+20./60.+49.24/3600)*15.,-42.-35./60.-52.5/3600.],
			 	[(13.+20./60.+27.01/3600)*15.,-42.-47./60.-29.0/3600.],
			 	[(13.+12./60.+42.78/3600)*15.,-42.-46./60.-52.2/3600.],
			 	[(13.+12./60.+10.10/3600)*15.,-42.-46./60.-48.8/3600.],
			 	[(13.+11./60.+45.06/3600)*15.,-42.-47./60.-26.5/3600.],
			 	[(13.+17./60.+48.98/3600)*15.,-42.-55./60.-45.8/3600.],
			 	[(13.+19./60.+01.98/3600)*15.,-43.-07./60.-49.5/3600.], #could be something special
			 	[(13.+18./60.+15.63/3600)*15.,-43.-06./60.-39.2/3600.],  #could be something special
			 	[(13.+26./60.+48.37/3600)*15.,-43.-41./60.-39.1/3600.],
			 	[(13.+24./60.+56.05/3600)*15.,-43.-40./60.-14.0/3600.],
			 	[(13.+27./60.+03.90/3600)*15.,-43.-19./60.-47.8/3600.],
			 	[(13.+24./60.+31.00/3600)*15.,-43.-28./60.-13.6/3600.],
			 	[(13.+22./60.+52.70/3600)*15.,-43.-20./60.-02.6/3600.],
			 	[(13.+23./60.+10.87/3600)*15.,-43.-01./60.-51.5/3600.],
			 	[(13.+26./60.+15.87/3600)*15.,-42.-47./60.-06.9/3600.],
			 	[(13.+26./60.+34.87/3600)*15.,-42.-43./60.-08.9/3600.],
			 	[(13.+24./60.+36.95/3600)*15.,-42.-34./60.-02.7/3600.],
			 	[(13.+23./60.+33.10/3600)*15.,-42.-15./60.-11.8/3600.],
			 	[(13.+27./60.+44.83/3600)*15.,-42.-06./60.-13.7/3600.]]
#			 	[13.27:43.23 -41:28:54.6],
#			 	[]
			 	
probes = [[cena_ra-0.1,cena_dec+0.1],[198.0,-45.1],[205.0,-41.0]]
probes = np.array(probes)

#gc_ra, gc_dec = np.loadtxt('wood2007_gc_coords.txt',unpack=True,comments='#')
#num, gc_ra, gc_dec = np.loadtxt('/Users/matt/Projects/081.D-0651A/data/CenA_allGC_degcoords.txt',unpack=True,comments='#')

pres = 0.339 #arcsec per pixel

gap_sht = 4.9  #arcmin
gap_lng = 10.4 #arcmin
chip_x = 2048. #pixels
chip_y = 2048. #pixels

gap_sht = gap_sht/60. #degrees
gap_lng = gap_lng/60. #degrees
chip_x = chip_x*pres/3600. #degrees
chip_y = chip_y*pres/3600. #degrees

fov_ra = 1.292 #degrees
fov_dec = 1.017 #degrees

hits_ra = []
hits_dec = []

fileout = open('pointing_coords.txt','w')


fig, ax = plt.subplots(figsize=(12,10))

n = 20. #Number of pointings in the spiral
d = 1.  #Diameter of spiral

#Get the coordinates for all the tiles
tiles = get_tiles()

#Cycle through 5 sets of 2 tiles and determine the spiral pointings for each set of two tiles
spiral1 = get_spiral(tiles[0])
spiral2 = get_spiral(tiles[1])

# print spiral1
# print spiral2
plt.figure()
for ii in range(len(spiral1)):
	plt.plot([spiral1[ii][0]],[spiral1[ii][1]],'bo')
	plt.plot([spiral2[ii][0]],[spiral2[ii][1]],'ro')

plt.xlim(cena_ra+2.,cena_ra-2.)
plt.ylim(cena_dec-2.,cena_dec+2.)
# plt.show()

#For each spiral pointing, find coordinates of each of the 6 pawprints
#e.g. Tile1, spiral pointing1, pawprint1 -> Tile2, spiral pointing1, pawprint1 -> T1, sp1, paw2 -> T2, sp1, paw2 -> ... -> T1, sp1,paw6 -> T2, sp1, paw6 -> T1, sp2, paw1 -> T2, sp2, paw1 -> ... -> T9, spN, paw6 -> t10, spN, paw6
coords_final = get_paw(spiral1,spiral2)

# print coords_final
# print len(coords_final)
plt.figure()
for ii in range(12):#len(coords_final)):
	plt.plot([coords_final[ii][0]],[coords_final[ii][1]],'ko')

plt.xlim(cena_ra+2.,cena_ra-2.)
plt.ylim(cena_dec-2.,cena_dec+2.)

# p1 = draw_chips(coords_final[0][0],coords_final[0][1],chip_x,gap_lng)
collection = []
for ii in range(12):#len(coords_final)):
	collection.append(draw_chips(coords_final[ii][0],coords_final[ii][1],chip_x,gap_lng))
	ax.plot([coords_final[ii][0]],[coords_final[ii][1]],'kx',ms=3)

for patches in collection:
	ax.add_collection(patches)

for ii in range(len(spiral1)):
	ax.plot([spiral1[ii][0]],[spiral1[ii][1]],'ko')
	ax.plot([spiral2[ii][0]],[spiral2[ii][1]],'ko')
	ax.annotate('%i' % (ii+1),(spiral1[ii][0],spiral1[ii][1]))
	ax.annotate('%i' % (ii+1),(spiral2[ii][0],spiral2[ii][1]))


#ax.plot(points_x,points_y,'ko')
ax.set_xlim(cena_ra+2.,cena_ra-2.)
ax.set_ylim(cena_dec-2.,cena_dec+2.)

plt.show()



exit()

#dra = 60./3600.
#ddec = 75./3600.

dra = 41./3600.
ddec = 54.0/3600.

#print dra*3600., ddec*3600.
#exit()

pointings_all_dec = []
for ii in range(len(pointings)):
	points = [[pointings[ii][0],pointings[ii][1]],
			  [pointings[ii][0]-1.*dra,pointings[ii][1]-1.*ddec],
			  [pointings[ii][0]+dra,pointings[ii][1]-1.*ddec],
			  [pointings[ii][0]+dra,pointings[ii][1]+ddec],
			  [pointings[ii][0]-1.*dra,pointings[ii][1]+ddec]]
	print >> fileout, "p%i" % (ii+1)
	for jj in range(len(points)):
		print >> fileout, points[jj][0], points[jj][1]

fileout.close()

for ii in range(len(pointings)):
	point_ra = pointings[ii][0]
	point_dec = pointings[ii][1]	

	print point_ra, point_dec

	points_x, points_y, collection = draw_footprint(point_ra,point_dec,chip_x,gap_sht)

	for patches in collection:
		ax.add_collection(patches)

	ax.plot(points_x,points_y,'ko')
	ax.annotate('%i' % (ii+1),(pointings[ii][0]-0.1,pointings[ii][1]+0.1))

print "Total hits = %i" % len(hits_ra)
hits_ra1 = []
hits_ra0 = []
hits_dec0 = []
hits_dec1 = []
for ii in range(len(hits_ra)):
	hits_ra0.append(hits_ra[ii][0])
	hits_ra1.append(hits_ra[ii][1])
	hits_dec0.append(hits_dec[ii][0])
	hits_dec1.append(hits_dec[ii][1])

pcheck_ra = np.linspace(np.min(hits_ra0),np.max(hits_ra1),1e4)
pcheck_dec = np.linspace(np.min(hits_dec0),np.max(hits_dec1),1e4)

counts = []

if sys.argv[2] == 'probe':
	print "Probing coverage..."
	for ii in range(len(probes)):
		print "Probe %i, RA...." % (ii+1)
		cnt_ra = []	
		for ra in pcheck_ra:
			mask = np.where(hits_dec0 <= probes[ii][1],1,0)
			chk_ra0 = np.compress(mask,hits_ra0)
			chk_ra1 = np.compress(mask,hits_ra1)
			chk_dec0 = np.compress(mask,hits_dec0)
			chk_dec1 = np.compress(mask,hits_dec1)
			mask = np.where(chk_dec1 >= probes[ii][1],1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_ra0 <= ra,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_ra1 >= ra,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
		
			cnt_ra.append(len(chk_ra0))
	
		print "Probe %i, Dec...." % (ii+1)
		cnt_dec = []	
		for dec in pcheck_dec:
			mask = np.where(hits_ra0 <= probes[ii][0],1,0)
			chk_ra0 = np.compress(mask,hits_ra0)
			chk_ra1 = np.compress(mask,hits_ra1)
			chk_dec0 = np.compress(mask,hits_dec0)
			chk_dec1 = np.compress(mask,hits_dec1)
			mask = np.where(chk_ra1 >= probes[ii][0],1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_dec0 <= dec,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
			mask = np.where(chk_dec1 >= dec,1,0)
			chk_ra0 = np.compress(mask,chk_ra0)
			chk_ra1 = np.compress(mask,chk_ra1)
			chk_dec0 = np.compress(mask,chk_dec0)
			chk_dec1 = np.compress(mask,chk_dec1)
		
			cnt_dec.append(len(chk_ra0))
	
		counts.append([cnt_ra,cnt_dec])
	
	counts = np.array(counts)

theta = np.linspace(0,2*np.pi,1e5)
r = np.arctan(3.e5/3.8e6)*180./np.pi
y = cena_dec+r*np.sin(theta)
x = cena_ra-r*np.cos(theta)/np.cos(y*np.pi/180.)
print r, cena_ra, cena_dec
print np.min(x), np.max(x)
print np.min(y), np.max(y)

ax.plot(x,y,'r--')

wcen, = ax.plot([wcen_ra],[wcen_dec],'o',color='blue',ms=20)

#for ii in range(len(probes)):
#	ax.plot([xlim[0],xlim[1]],[probes[ii][1],probes[ii][1]],'k--')
#	ax.plot([probes[ii][0],probes[ii][0]],[ylim[0],ylim[1]],'k--')
#	ax.plot([probes[ii][0]],[probes[ii][1]],'kx',ms=20)

for ii in range(len(gal_coords)):
	dwarf, = ax.plot([gal_coords[ii][0]],[gal_coords[ii][1]],'o',color='darkorange',ms=15,markeredgecolor=None)#,alpha=0.75)

for ii in range(len(crno_coords)):
	crno, = ax.plot([crno_coords[ii][0]],[crno_coords[ii][1]],'p',color='maroon',ms=12,markeredgecolor=None)#alpha=0.75)
gcs, = ax.plot(gc_ra,gc_dec,'k,')

cena, = ax.plot([cena_ra],[cena_dec],'*',color='green',ms=25)

ax.legend((cena,gcs,dwarf,crno,wcen),(r'NGC$\,$5128',r'Confirmed GCs $$',r'Known Dwarfs $$',r'Dwarf Candidates $$',r'$\omega$ Centauri'),loc=3,numpoints=1)

ax.set_xlim(xlim[0],xlim[1])
ax.set_ylim(ylim[0],ylim[1])
# ax.set_xlim(208.8,195.0)
# ax.set_ylim(ylim[0],ylim[1])
ax.tick_params(axis='both',labelsize=18)
ax.set_xlabel(r'$\alpha$ [deg]',fontsize=20)
ax.set_ylabel(r'$\delta$ [deg]',fontsize=20)

plt.savefig('dithering_pattern_4jobs.jpg')
#plt.savefig('dithering_pattern_4poster.pdf')

if sys.argv[2] == 'probe':
	plt.figure(figsize=(10,10))
	cnt = 0
	for ii in range(len(probes)):
		cnt = cnt + 1
		plt.subplot('%i%i%i' % (len(probes),2,cnt))
		plt.plot(pcheck_ra,counts[ii][0]*dit,'b')
		plt.xlim(xlim[0],xlim[1])
		plt.xlabel(r'RA [deg]')
		plt.ylabel(r'Exp. Time [s]')
	
		cnt = cnt + 1
		plt.subplot('%i%i%i' % (len(probes),2,cnt))
		plt.plot(pcheck_dec,counts[ii][1]*dit,'b')
		plt.xlim(ylim[0],ylim[1])
		plt.xlabel(r'Dec [deg]')
		plt.ylabel(r'Exp. Time [s]')
	
	plt.savefig('decam_probes.pdf')

plt.show()
	