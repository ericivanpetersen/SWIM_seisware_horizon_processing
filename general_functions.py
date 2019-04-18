
import numpy as np
from osgeo import gdal
import csv
from mars_projections import *
import matplotlib
import matplotlib.pyplot as plt

def write_csv_file(outfile, data, header=None):
	"""Simple function to write csv file based on 
	header information and data matrix	
	
	:param str header: header indicating fields 
		contained withing csv file
	:param float data: array filled with data to 
		write to csv file
	"""
	
	with open(outfile, 'wb') as fout:
		writer = csv.writer(fout)
		if header:    writer.writerow(header)
		writer.writerows(data)

def extract_raster_profile(x, y, rasterfile):
	"""Function that extracts values from the raster
	at given X & Y values; assumes single band for
	the raster. X/Y values must be in same projection
	as the raster.

	:param str rasterfile: gdal raster file to read
	:param float x: x value at which to extract 
		raster value; can be an array
	:param float y: y value at which to extract
		raster value; can be an array
	:output float z: z value extracted from raster
	"""

	layer = gdal.Open(rasterfile)
	gt = layer.GetGeoTransform() # get coord to raster conversion
	# convert coordinates to raster:
	rasterx = np.rint( (x - gt[0]) / gt[1] )
	rastery = np.rint( (y - gt[3]) / gt[5] )
	# find values and return:
	NN = len(rasterx)
	z = np.empty((NN,1), dtype='float')
	for nn in range(NN):
		z[nn] = layer.GetRasterBand(1).ReadAsArray(rasterx[nn], rastery[nn], 1, 1)
	return z

def extract_MOLA_profile(lon, lat, rasterfile='/Users/eric/Documents/orig/supl/MOLA/DEM_global_mola128ppd_merged_mola64ppd/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif'):
	"""Extract topographic profile of MOLA elevations
	based on input Latitude & Longitude values.

	:param float lat: latitude in Mars 2000
	:param float lon: longitude in Mars 2000
	:param str rasterfile: path to MOLA geotif file; default
		is the location on Eric's computer
	:output float z: elevation value in meters
	"""

	# Convert lat/lon values to Mars Equirectangular:
	x, y = transform_XY(lon, lat, mars_2000, mars_equidistant_cylindrical)
	
	# Extract elevation values:
	z = extract_raster_profile(x, y, rasterfile)

	return z
	
def depth_correct_radar(TWT_surf, TWT_sub, epsilon):
	"""Depth-corrects radar data in twt (ns)
	using input dielectric constant epsilon

	:param float TWT_surf: two-way travel time to surface reflector
	:param float TWT_sub: two-way travel time to subsurface reflector
	:param float epsilon: assumed dielectric constant epsilon

	:output float depth: depth in meters to subsurface reflector
	"""
	
	depth = 0.3 / np.sqrt(epsilon) * (TWT_sub - TWT_surf) / 2

	return depth

def estimate_epsilon(TWT_surf, TWT_sub, depth):
	"""Estimates dielectric constant based on input TWT_surf,
	TWT_sub, and target depth

	:param float TWT_surf: two way travel time to surface reflector
	:param float TWT_sub: two-way travel time to subsurface reflector
	:param float depth: constrained depth to reflector/thickness of 	
		deposit	

	:output float epsilon: estimated dielectric constant
	"""

	epsilon = ( 0.3 * (TWT_sub - TWT_surf) / (2 * depth) ) ** 2

	return epsilon

def estimate_epsilon_targetTWT(TWT_surf, TWT_sub, TWT_targ):
	"""Estimates dielectric constant based on input 'target TWT'

	:param float TWT_surf: TWT to surface refl
	:param float TWT_sub: TWT to subsurface refl
	:param flaot TWT_targ: target TWT for depth converstion

	:output float epsilon: estimated dielectric constant
	"""

	depth = 0.3 * ( TWT_targ - TWT_surf ) / 2 # Synthetic depth/thickness based on target TWT
	epsilon = estimate_epsilon( TWT_surf, TWT_sub, depth)

	return epsilon

def estimate_epsilon_MOLA_minima_extrapolation(orbit, trace, lat, lon, TWT_surf, TWT_sub, MOLA_file, outpath):
	"""Estimates the dielectric constant for subsurface 
	reflectors on an individual SHARAD line using a method
	developed by Ali Bramson: MOLA elevation is extracted
	from MOLA geotiff to produce a topographic profile for
	the SHARAD line; the function then asks for user input on
	the range of latitudes for which the local minima is 
	expected to coincide with the intersection of the mapped
	subsurface reflector with the surface. Two minima, on
	either side of the mapped reflector, are selected and 
	then used to create a linear extrapolation to force the
	subsurface reflector to, estimating the resultant 
	dielectric constant.
	NOTE: assumes lon, lat, TWT_surf, and TWT_sub are the 
		same length; no-data points are =-9999.99

	:param float lon: longitude for each SHARAD trace
	:param float lat: latitude for each SHARAD trace
	:param float TWT_surf: two-way travel time for the 
		surface reflector
	:param float TWT_sub: two-way travel time for the 
		subsurface reflector
	:param float orbit: orbit number, for labelling the plot
	:param float trace: trace numbers, for indexing
	:param str MOLA_file: filepath to MOLA geotiff for 
		obtaining elevation profile
	:param str outpath: outpath to save plot
	:output float epsilon: estimated dielectric constant
		using the method
	:output float z_sub: the assumed z_sub elevations
		derived from the MOLA minima method used to 
		estimate the dielectric constant
	"""

	TWT_sub[(TWT_sub == -9999.99)] = np.nan # set NaNs for TWT_sub

	non_result = np.full(np.size(TWT_sub), np.nan) # set non-result
	
	# produce the MOLA profile:
	z_surf = extract_MOLA_profile(lon, lat, MOLA_file)
	# produce non-depth-corrected z_sub:
	z_sub_eps1 = z_surf.transpose() - depth_correct_radar(TWT_surf, TWT_sub, 1)
	z_surf = z_surf.squeeze()
	z_sub_eps1 = z_sub_eps1.squeeze()
	z_sub = np.full(np.size(z_surf), np.nan)
	epsilon = np.full(np.size(z_surf), np.nan)
	print
	print('Dielectric Estimation for SHARAD Line {}').format(orbit)
	print('Examine graph and determine # of desired line segments')
	print('        to define between observed minima in the MOLA profile,')
	print('	       as well as the range in latitude for each minima.')
	# plot the profile:
	plt.figure(1)
	plt.subplot(211)
	plt.title(str(orbit)+': Dielectric Estimation')
	plt.plot( lat, z_surf, lat, z_sub_eps1 )
	plt.ylabel('MOLA Elevation (m)')
	plt.legend(['Surface', 'Subsurface, eps=1'])
	plt.subplot(212)
	plt.plot( lat, -TWT_surf, lat, -TWT_sub)
	plt.ylabel('TWT (ns)')
	plt.legend(['Surface', 'Subsurface'])
	plt.xlabel('Latitude')
	plt.show()

	# ask for inputs on range to search for minima:
	num_seg = input("Number of Segments Between Minima? (if profile insufficient to estimate minima, enter '0') ")
	if num_seg == 0:
		return non_result, non_result, []
	# initialize segment summary array, which includes statistical
	#	info, center locations, and locations of elevation minima
	seg_sum = np.full((num_seg,11), np.nan)

	nn = 0
	while nn < num_seg:
		print
		print('Input minima ranges for line segment {}').format(nn+1)
		lat_min_1a = input("Start point to search for first minima: ")
		# if entered "0", exit the program, return nans:
		if lat_min_1a == 0 :
			return non_result, non_result, []
		lat_min_1b = input("  End point to search for first minima: ")
		lat_min_2a = input("Start point to search for second minima: ")
		lat_min_2b = input("  End point to search for second minima: ")
		# find desired local minima in elevation for those latitude ranges:
		minrange1 = np.where( (lat > lat_min_1a) & (lat < lat_min_1b) )
		minrange2 = np.where( (lat > lat_min_2a) & (lat < lat_min_2b) )
		zmin1 = np.amin(z_surf[minrange1])
		zmin2 = np.amin(z_surf[minrange2])
		zmin1_ind = np.where( (z_surf == zmin1) & (lat>lat_min_1a) & (lat<lat_min_1b))
		zmin2_ind = np.where( (z_surf == zmin2) & (lat>lat_min_2a) & (lat<lat_min_2b))
		zmin1_ind = zmin1_ind[0][-1]
		zmin2_ind = zmin2_ind[0][0]
		if zmin1_ind > zmin2_ind:
			seg_ind = np.arange(zmin2_ind, zmin1_ind+1, 1)
		else:
			seg_ind = np.arange(zmin1_ind, zmin2_ind+1, 1)
		lat1 = lat[zmin1_ind]
		lat2 = lat[zmin2_ind]
		# Find the center of the area where reflectors are mapped:
		trace_seg = trace[seg_ind]
		seg_refl_ind = np.where(z_sub[seg_ind] != np.nan)
		trace1 = np.amin(trace_seg[seg_refl_ind])
		trace2 = np.amax(trace_seg[seg_refl_ind])
		center_trace = np.around( (trace2-trace1)/2 ) + trace1
		center_lon = lon[np.where(trace == center_trace)]
		center_lat = lat[np.where(trace == center_trace)]
		# calculate target depths from linear extrapolation between minima:
		slope = (zmin2 - zmin1) / (lat2 - lat1)
		z_sub[seg_ind] = zmin1 + slope * (lat[seg_ind] - lat1)
		# calculate estimated dielectric constant epsilon:
		epsilon[seg_ind] = estimate_epsilon(TWT_surf[seg_ind], TWT_sub[seg_ind], (z_surf[seg_ind]-z_sub[seg_ind]))
		# Plot the profile again, with results:
		print('Median Epsilon = {}').format(np.nanmedian(epsilon[seg_ind]))
		fig = plt.figure(2)
		plt.subplot(211)
		plt.title(str(orbit)+': Real Dielectric Estimation')
		plt.plot( lat, z_surf, lat, z_sub_eps1, lat, z_sub)
		plt.legend(['Surface', 'Subsurface, eps=1', 'Target Base'])
		plt.ylabel('MOLA Elevation (m)')
		plt.axis([np.amin(lat), np.amax(lat), np.nanmin(z_sub_eps1)-50, np.nanmax(z_surf)+50])
		plt.subplot(212)
		plt.plot( lat, epsilon, '*')
		plt.ylabel('Dielectric Constant')
		plt.axis([np.amin(lat), np.amax(lat), 1, 12])
		plt.xlabel('Latitude')
		plt.show()

		result = raw_input("Satisfied with segment? (y=yes, n=no, r=retry, ra=retry all) ")
		if result == 'y':
			# Save segment summary results:
			seg_sum[nn,0] = orbit
			seg_sum[nn,1] = center_lat
			seg_sum[nn,2] = center_lon
			seg_sum[nn,3] = np.nanmedian(epsilon[seg_ind])
			seg_sum[nn,4] = np.nanmean(epsilon[seg_ind])
			seg_sum[nn,5] = np.nanstd(epsilon[seg_ind])
			seg_sum[nn,6] = np.subtract(*np.nanpercentile(epsilon[seg_ind], [75, 25]))
			seg_sum[nn,8] = lon[zmin1_ind]
			seg_sum[nn,7] = lat1
			seg_sum[nn,10] = lon[zmin2_ind]
			seg_sum[nn,9] = lat2
			# Proceed
			nn = nn+1
			# Save figure & Return
			if nn == num_seg: 
				fig.savefig(outpath+str(orbit)+'.pdf')
				return epsilon, z_sub, seg_sum
		if result == 'n':
			epsilon[seg_ind] = np.full(np.size(seg_ind), np.nan)
			z_sub[seg_ind] = np.full(np.size(seg_ind), np.nan)
			nn = nn+1
			if nn == num_seg:
				return epsilon, z_sub, seg_sum
		if result == 'r':
			pass
		if result == 'ra':
			nn = 0

def convert_epsilon_to_conf(epsilon):
	"""Calculates SWIM consistency value (con) from epsilon
	based on our transfer function where 1 is eps=3, 
	0 is eps=5, & -1 is eps=7

	:param float epsilon: dielectric constant
	
	:output float con: consistency value between -1 & 1
	"""

	conf = 2.5 - 0.5 * epsilon
	
	# Restrict to values between -1 & 1:
	conf[np.argwhere(conf > 1)] = 1
	conf[np.argwhere(conf < -1)] = -1

	# Extreme epsilon values resolve to con = 0:
	conf[np.argwhere(epsilon < 1)] = 0
	conf[np.argwhere(epsilon > 12)] = 0 	

	return conf


def convert_xyz_to_raster(x, y, z, res, method, outfile):
	"""Function to produce raster from point data
	presented with X, Y, Z Coordinates

	NOTE: This function still has some issues,
		particularly with outputting large rasters
		(which have to be cut into smaller chunks).

	:param float x: x coordinates
	:param float y: y coordinates
	:param float z: z coordinates or other variable
		that you want to display in the raster
	:param float res: resolution of the output 
		raster in relevant projection units
	:param str outfile: file to save raster.
		NOTE: if size of raster is large, outfile
		will be split into numbered raster files
	"""

	bean = 500 # Minimum size (square) of output raster before cutting it into smaller chunks 
		# This is to resolve an issue with the slowness of gdal.Grid

	# Determine size of output grid based on desired resolution
	w1 = np.ceil( np.ptp(x) / res )
	h1 = np.ceil( np.ptp(y) / res)
	# Define the algorithm for extrapolating data:
	alg = method+':radius1='+str(res/2)+':radius2='+str(res/2)+':nodata=-9999'

	# As long as output raster size is less than bean x bean, output one raster:
	if w1 * h1 < bean**2:
	
		xyz = np.column_stack((x, y, z))
		# Sort by x, then by y:
		xyz = xyz[xyz[:,0].argsort()]
		xyz = xyz[xyz[:,1].argsort(kind='mergesort')]
		
		write_csv_file('grid.csv',xyz) # Write temporary xyz file 
		#Convert to Grid
		gdal.Grid(outfile,'grid.vrt', width=w1, height=h1, algorithm=alg)
		os.remove('./grid.csv') # Remove temporary xyz file
		
		print 
		print('Produced single map with filename {}'.format(outfile))
		print
	# If size larger than bean x bean, split large grid into subgrids for faster processing:
	else:
		#merge_command = ['', '-o', outfile] #attempting merging rasters; this never worked out
		# Cutting into squares size bean x bean:
		Wind = int(np.ceil(w1/bean))
		Hind = int(np.ceil(h1/bean))
		nn = 0
		for W in range(Wind):
			x1 = W*bean*res + np.min(x)
			x2 = (W+1)*bean*res + np.min(x)
			xind = (x >= x1) & (x < x2)
			if W < Wind:
				wt = bean
			else:
				wt = w1 % bean
			for H in range(Hind):
				y1 = (H)*bean*res + np.min(y)
				y2 = (H+1)*bean*res + np.min(y)
				yind = (y >= y1) & (y < y2)
				ind = xind & yind
	
				# Check for data within current square:
				if any(ind):
					xyz_temp = np.column_stack((x[ind], y[ind], z[ind]))
					xyz_temp = xyz_temp[xyz_temp[:,0].argsort()]
					xyz = xyz_temp[xyz_temp[:,1].argsort(kind='mergesort')]
					if H < Hind:
						ht = bean
					else:
						ht = h1 % bean

					out = outfile[:-4] + str(nn+1) + '.tif' # Numbered outfile
					write_csv_file('grid.csv',xyz_temp)
					gdal.Grid(out, 'grid.vrt', width=wt, height=ht, algorithm=alg)
					os.remove('./grid.csv') #Remove temporary xyz file
					#merge_command.append(out)  
					nn = nn+1
		print
		print 'Dataset too large for single grid production'
		print('Produced {} files with the following names: '.format(nn))
		for n in range(nn-1):
			print('{}'.format(outfile[:-4] + str(n+1) + '.tif'))
		print

		# Abandoned Merge Raster Code:
		#num_files = nn
		#print(merge_command)
		#gdal_merge.main(merge_command)		

		#for NN in range(num_files):
		#	temp = './temp'+str(NN)+'.tif'
		#	os.remove(temp)

if __name__ == '__main__':
	lon = 20 + np.arange(11)
	lat = 40 + np.arange(11)
	#d = np.sqrt( x**2 + y**2)
	z = extract_MOLA_profile(lat, lon)
	print z
