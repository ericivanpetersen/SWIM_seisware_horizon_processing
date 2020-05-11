
import numpy as np
from osgeo import gdal
import csv
from mars_projections import *
import matplotlib.pyplot as plt

def write_csv_file(outfile, data, header=None):
	"""Simple function to write csv file based on 
	header information and data matrix	
	
	:param str header: header indicating fields 
		contained withing csv file
	:param float data: array filled with data to 
		write to csv file
	"""
	
	with open(outfile, 'w') as fout:
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

def extract_MOLA_profile(lon, lat, rasterfile='../mola_data/dem/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif'):
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
	# Squeeze out unnecessary dimensions:
	z_surf = z_surf.squeeze()
	z_sub_eps1 = z_sub_eps1.squeeze()
	# Create arrays to fill with results:
	z_sub = np.full(np.size(z_surf), np.nan)
	epsilon = np.full(np.size(z_surf), np.nan)

	# Begin Initial Assessment Step:
	print()
	print('Dielectric Estimation for SHARAD Line {}'.format(orbit))
	print('Examine graph and determine # of desired line segments')
	print('        to define between observed minima in the MOLA profile,')
	print('	       as well as the range in latitude for each minima.')
	# plot the profile:
	plt.figure(1)
	plt.subplot(211)
	plt.title(str(orbit)+': Dielectric Estimation')
	plt.plot( trace, z_surf, trace, z_sub_eps1 )
	plt.ylabel('MOLA Elevation (m)')
	plt.legend(['Surface', 'Subsurface, eps=1'])
	plt.subplot(212)
	plt.plot( trace, -TWT_surf, trace, -TWT_sub)
	plt.ylabel('TWT (ns)')
	plt.legend(['Surface', 'Subsurface'])
	plt.xlabel('Trace')
	plt.show()

	# ask for inputs on range to search for minima:
	num_seg = int(input("Number of Segments Between Minima? (if profile insufficient to estimate minima, enter '0') "))
	if num_seg == 0:
		return non_result, non_result, []
	# initialize segment summary array, which includes statistical
	#	info, center locations, and locations of elevation minima
	seg_sum = np.full((num_seg,11), np.nan)

	nn = 0
	while nn < num_seg:

		# Plot the profile & request user input:
		print()
		print("A plot will be shown of picks from the radargram")
		print("  For a segment of subsurface picks, select four points:")
		print("   --Two for each minima, selecting the region in which you wish to seek the minima")
		plt.figure(1)
		plt.title(str(orbit)+': Dielectric Estimation')
		plt.plot( trace, z_surf, trace, z_sub_eps1 )
		plt.ylabel('MOLA Elevation (m)')
		plt.legend(['Surface', 'Subsurface, eps=1'])
		plt.xlabel('Trace')
		tlist = plt.ginput(n=4) # Asks for input on plot:
		plt.close(1)
		# Extract Trace Values for Selected Points:
		ttt = [tlist[0][0], tlist[1][0], tlist[2][0], tlist[3][0]]
		ttt.sort()
		t_min_1a = ttt[0]
		t_min_1b = ttt[1]
		t_min_2a = ttt[2]
		t_min_2b = ttt[3]

		# find desired local minima in elevation for those latitude ranges:
		minrange1 = np.where( (trace > t_min_1a) & (trace < t_min_1b) )
		minrange2 = np.where( (trace > t_min_2a) & (trace < t_min_2b) )
		zmin1 = np.amin(z_surf[minrange1])
		zmin2 = np.amin(z_surf[minrange2])
		zmin1_ind = np.where( (z_surf == zmin1) & (trace>t_min_1a) & (trace<t_min_1b))
		zmin2_ind = np.where( (z_surf == zmin2) & (trace>t_min_2a) & (trace<t_min_2b))
		zmin1_ind = zmin1_ind[0][-1]
		zmin2_ind = zmin2_ind[0][0]
		if zmin1_ind > zmin2_ind:
			seg_ind = np.arange(zmin2_ind, zmin1_ind+1, 1)
		else:
			seg_ind = np.arange(zmin1_ind, zmin2_ind+1, 1)

		# Find the center of the area where reflectors are mapped:
		trace_seg = trace[seg_ind]
		seg_refl_ind = np.invert(np.isnan(TWT_sub[seg_ind]))
		trace1 = np.amin(trace_seg[seg_refl_ind])
		trace2 = np.amax(trace_seg[seg_refl_ind])
		center_trace = np.around( (trace2-trace1)/2 ) + trace1
		center_ind = np.where( np.absolute( trace - center_trace ) == np.amin(np.absolute( trace - center_trace))) #Closest to center, in case that trace has been skipped. 
		center_lon = lon[center_ind]
		center_lat = lat[center_ind]
		center_lon = center_lon[0]
		center_lat = center_lat[0]

		# calculate target depths from linear extrapolation between minima:
		slope = (zmin2 - zmin1) / (trace2 - trace1)
		z_sub[seg_ind] = zmin1 + slope * (trace[seg_ind] - trace1)
		# calculate estimated dielectric constant epsilon:
		epsilon[seg_ind] = estimate_epsilon(TWT_surf[seg_ind], TWT_sub[seg_ind], (z_surf[seg_ind]-z_sub[seg_ind]))
		# Plot the profile again, with results:
		print('Median Epsilon = {}'.format(np.nanmedian(epsilon[seg_ind])))
		fig = plt.figure(2)
		plt.subplot(211)
		plt.title(str(orbit)+': Real Dielectric Estimation')
		plt.plot( trace, z_surf, trace, z_sub_eps1, trace, z_sub)
		plt.legend(['Surface', 'Subsurface, eps=1', 'Target Base'])
		plt.ylabel('MOLA Elevation (m)')
		plt.axis([np.amin(lat), np.amax(lat), np.nanmin(z_sub_eps1)-50, np.nanmax(z_surf)+50])
		plt.subplot(212)
		plt.plot( trace, epsilon, '*')
		plt.ylabel('Dielectric Constant')
		plt.axis([np.amin(lat), np.amax(lat), 1, 12])
		plt.xlabel('Trace')
		plt.show()

		result = input("Satisfied with segment? (y=yes, n=no, r=retry, ra=retry all) ")
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
			seg_sum[nn,7] = lat[zmin1_ind]
			seg_sum[nn,10] = lon[zmin2_ind]
			seg_sum[nn,9] = lat[zmin2_ind]
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

def estimate_epsilon_MOLA_polyfit(orbit, trace, lat, lon, TWT_surf, TWT_sub, MOLA_file, outpath, pn=1):
	"""Estimates the dielectric constant for subsurface 
	reflectors on an individual SHARAD line by fitting a 
	polynomial (n=1, linear by default) to a selected section
	of surface reflector and extending that function to the 
	
	The function then asks for user input on
	the range of latitudes between which to fit an
	Nth-order polynomial on MOLA elevation, extrapolating
	beneath the feature of interest (LDA, etc) to force the
	observed subsurface reflector to, estimating the resultant 
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
	:param int pn: order of the polynomial fit.
	:output float epsilon: estimated dielectric constant
		using the method
	:output float z_sub: the assumed z_sub elevations
		derived from the MOLA minima method used to 
		estimate the dielectric constant
	"""

	TWT_sub[(TWT_sub == -9999.99)] = np.nan # set NaNs for TWT_sub
	non_result = np.full(np.size(TWT_sub), np.nan) # set non-result
	
	#Produce surface profiles:
	z_surf = -TWT_surf * 0.3/2 # SHARAD-only surface elevation profile
	z_surf_MOLA = extract_MOLA_profile(lon, lat, MOLA_file) # MOLA profile extracted from nadir points
	# Produce non-depth-corrected z_sub (epsilon = 1) profiles:
	z_sub_eps1 = -TWT_sub * 0.3/2 
	z_sub_eps1_MOLA = z_surf_MOLA.transpose() - depth_correct_radar(TWT_surf, TWT_sub, 1)
	# Squeeze out unnecesssary dimensions:
	z_surf = z_surf.squeeze()
	z_sub_eps1 = z_sub_eps1.squeeze()
	z_surf_MOLA = z_surf_MOLA.squeeze()
	z_sub_eps1_MOLA = z_sub_eps1_MOLA.squeeze()

	# Create arrays to fill with results:
	z_sub = np.full(np.size(z_surf), np.nan)
	z_sub_MOLA = np.full(np.size(z_surf), np.nan)
	z_sub_flat = np.full(np.size(z_surf), np.nan)
	epsilon = np.full(np.size(z_surf), np.nan)
	epsilon_MOLA = np.full(np.size(z_surf), np.nan)
	eps_flat = np.full(np.size(z_surf), np.nan)
	# Initial assessment step:
	print()
	print('Dielectric Estimation for SHARAD Line {}'.format(orbit))
	print('Examine graph and determine # of desired line segments')
	print('        in the subsurface reflector for which you''d like to')
	print('	       define an extrapolated depth-correction surface.')
	# Plot the profile:
	plt.figure(1)
	plt.subplot(211)
	plt.title(str(orbit)+': Dielectric Estimation')
	plt.plot( trace, z_surf_MOLA, trace, z_sub_eps1_MOLA )
	plt.ylabel('Elevation (m)')
	plt.legend(['Surface', 'Subsurface, eps=1'])
	plt.subplot(212)
	plt.plot( trace, -TWT_surf, trace, -TWT_sub)
	plt.ylabel('TWT (ns)')
	plt.legend(['Surface', 'Subsurface'])
	plt.xlabel('Trace')
	plt.show()

	# Ask for inputs on number of individual segments to estimate dielectric for:
	num_seg = int(input("Number of Individual Segments to Split Reflector into (if profile insufficient to estimate epsilon, enter '0') "))
	if num_seg == 0:
		return non_result, non_result, [], non_result
	# initialize segment summary array, which includes statistical
	#	info, center locations, and locations of elevation minima
	seg_sum = np.full((num_seg,11), np.nan)

	nn = 0
	while nn < num_seg:

		# Plot the profile & request user input:
		print()
		print("A plot will be shown of picks from the radargram")
		print("  For a segment of subsurface picks, select three points:")
		print("   --First click on two points defining the start and end points")
		print("        for the region you want to fit the polynomial to")
		print("   --Then click on a third point defining the extent to which you")
		print("         wish to extend that fit for estimating dielectric constant")
		plt.figure(1)
		plt.title(str(orbit)+': Dielectric Estimation')
		plt.plot( trace, z_surf, trace, z_sub_eps1 )
		plt.ylabel('Elevation (m)')
		plt.legend(['Surface', 'Subsurface, eps=1'])
		tlist = plt.ginput(n=3) # Asks for input on plot:
		plt.close(1)
		# Extract Trace Values for Selected Points:
		t1 = tlist[0][0]
		t2 = tlist[1][0]
		t3 = tlist[2][0]
		
		# Organize input traces to produce fit & segment indices:
		t_fit_min = min(t1,t2)
		t_fit_max = max(t1,t2)
		t_ext_min = min(t1,t2,t3)
		t_ext_max = max(t1,t2,t3)
		if t_ext_max == t3:
			t_nearest = t_fit_max
		if t_ext_min == t3:
			t_nearest = t_fit_min
		z_nearest = z_surf[np.argmin(np.abs(trace-t_nearest))]
		fitind = np.where((trace >= t_fit_min) & (trace <= t_fit_max))
		seg_ind = np.where((trace >= t_ext_min) & (trace <= t_ext_max))

		# Perform polynomial fit and 
		pfit = np.polyfit(lat[fitind], z_surf[fitind], pn)
		z_sub[seg_ind] = np.polyval(pfit, lat[seg_ind])
		# Same for MOLA:
		pfit_MOLA = np.polyfit(lat[fitind], z_surf_MOLA[fitind], pn)
		z_sub_MOLA[seg_ind] = np.polyval(pfit_MOLA, lat[seg_ind])
		# calculate estimated dielectric constant epsilon:
		epsilon[seg_ind] = estimate_epsilon(TWT_surf[seg_ind], TWT_sub[seg_ind], (z_surf[seg_ind]-z_sub[seg_ind]))
		epsilon_MOLA[seg_ind] = estimate_epsilon(TWT_surf[seg_ind], TWT_sub[seg_ind], (z_surf_MOLA[seg_ind]-z_sub_MOLA[seg_ind]))
		
		# Calculate Dielectric Constant Required for Flat Surface (For Plot Comparison)
		z_sub_flat[seg_ind] = z_nearest
		eps_flat[seg_ind] = estimate_epsilon(TWT_surf[seg_ind], TWT_sub[seg_ind], (z_surf[seg_ind]-z_sub_flat[seg_ind]))
		# Find the center of the area where reflectors are mapped:
		trace_seg = trace[seg_ind]
		seg_refl_ind = np.invert(np.isnan(TWT_sub[seg_ind]))
		trace1 = np.amin(trace_seg[seg_refl_ind])
		trace2 = np.amax(trace_seg[seg_refl_ind])
		center_trace = np.around( (trace2-trace1)/2 ) + trace1
		center_ind = np.where( np.absolute( trace - center_trace ) == np.amin(np.absolute( trace - center_trace))) #Closest to center, in case that trace has been skipped. 
		center_lon = lon[center_ind]
		center_lat = lat[center_ind]
		center_lon = center_lon[0]
		center_lat = center_lat[0]
		
		# Plot the profile again, with results:
		print('Median Epsilon = {}'.format(np.nanmedian(epsilon[seg_ind])))
		print('Median Flat Epsilon = {}'.format(np.nanmedian(eps_flat[seg_ind])))
		fig = plt.figure(2)
		plt.subplot(211)
		plt.title(str(orbit)+': Real Dielectric Estimation')
		plt.plot( trace, z_surf, 'k')
		plt.plot(trace, z_sub_eps1, 'g')
		plt.plot(trace, z_sub_flat,'r')
		plt.plot(trace, z_sub, 'b')
		plt.legend(['Surface', 'Subsurface, eps=1', 'Flat Base', 'Sloped Base'])
		plt.ylabel('Elevation (m)')
		plt.axis([np.amin(lat), np.amax(lat), np.nanmin(z_sub_eps1)-50, np.nanmax(z_surf)+50])
		plt.subplot(212)
		plt.plot(trace, eps_flat, 'r*')
		plt.plot( trace, epsilon, 'b*')
		plt.ylabel('Dielectric Constant')
		plt.legend(['Flat Base', 'Sloped Base'])
		plt.axis([np.amin(trace), np.amax(trace), 1, 12])
		plt.xlabel('Trace')
		plt.show()

		result = input("Satisfied with segment? (y=yes, n=no, r=retry, ra=retry all) ")
		if result == 'y':
			# Save segment summary results:
			seg_sum[nn,0] = orbit
			seg_sum[nn,1] = center_lat
			seg_sum[nn,2] = center_lon
			seg_sum[nn,3] = np.nanmedian(epsilon[seg_ind])
			seg_sum[nn,4] = np.nanmean(epsilon[seg_ind])
			seg_sum[nn,5] = np.nanstd(epsilon[seg_ind])
			seg_sum[nn,6] = np.subtract(*np.nanpercentile(epsilon[seg_ind], [75, 25]))
			seg_sum[nn,7] = t1
			seg_sum[nn,8] = t2
			seg_sum[nn,9] = t3
			# Proceed
			nn = nn+1
			# Save figure & Return
			if nn == num_seg: 
				fig.savefig(outpath+str(orbit)+'.pdf')
				return epsilon, z_sub, seg_sum, epsilon_MOLA
		if result == 'n':
			epsilon[seg_ind] = np.full(np.size(seg_ind), np.nan)
			z_sub[seg_ind] = np.full(np.size(seg_ind), np.nan)
			nn = nn+1
			if nn == num_seg:
				return epsilon, z_sub, seg_sum, epsilon_MOLA
		if result == 'r':
			pass
		if result == 'ra':
			nn = 0

def estimate_epsilon_min_R2(trace, TWT_surf, TWT_sub, plot=False, eps_range=np.arange(1,8,0.01)):
	"""Estimates the best-fit dielectric constant
	to reduce the correlation between surface and
	subsurface reflectors.

	:param float trace: radargram trace number
	:param float TWT_surf: TWT for surface return
	:param float TWT_sub: TWT for subsurface return
	:param bool plot: indicates whether or not to plot
	:output float eps_best: best-fit epsilon value
	:output float depth: depth values for eps_best
	:output float r2: minimum r2 value
	"""

	EE = np.size(eps_range)
	r2_range = np.full(EE, 1.5)
	for ee in range(EE):
		tmp_sub_corr = TWT_surf + (TWT_sub-TWT_surf)/np.sqrt(eps_range[ee])
		tmp_depth = depth_correct_radar(TWT_surf, TWT_sub, eps_range[ee])
		r2_range[ee] = np.corrcoef(TWT_surf, tmp_sub_corr)[0,1]**2

	r2_min = np.min(r2_range)
	mind = np.where(r2_range == r2_min)[0][0]
	eps_best = eps_range[mind]
	sub_corr = TWT_surf + (TWT_sub-TWT_surf)/np.sqrt(eps_best)
	depth_best = depth_correct_radar(TWT_surf, TWT_sub, eps_best)
	sub3 = TWT_surf + (TWT_sub - TWT_surf)/np.sqrt(3)
	print("Min R2 = {0} at Eps = {1}".format(r2_min, eps_best))
	if plot==True:
		plt.rcParams.update({'font.size':16})

		plt.figure(0)
		#plt.subplot(1,2,1)
		plt.plot(eps_range, r2_range, 'k')
		plt.xlabel('$\epsilon$\'')
		plt.ylabel(r'R$^2$')
		plt.show()

		plt.figure(1)
		#plt.subplot(1,2,2)
		plt.plot(trace, -TWT_surf, trace, -TWT_sub, trace, -sub_corr, trace, -sub3)
		plt.xlabel('Trace')
		plt.ylabel('TWT (ns)')
		plt.legend(['Surface', 'Subsurface', 'Corrected', 'Eps=3'])
		plt.show()

		plt.figure(2)
		#plt.subplot(1,3,1)
		plt.plot(-TWT_surf, -TWT_sub, 'b*')
		plt.xlabel('Surface TWT (ns)')
		plt.ylabel('Subsurface TWT (ns)')
		#plt.subplot(1,3,2)
		plt.plot(-TWT_surf, -sub_corr, 'k*')
		#plt.subplot(1,3,3)
		plt.plot(-TWT_surf, -tmp_sub_corr, 'r*')
		#plt.xlabel('Overcorrected')
		plt.legend(['Uncorrected', 'Corrected', 'Overcorrected'])
		plt.show()

	return eps_best, depth_best, r2_min

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


if __name__ == '__main__':
	# Test MOLA profile extraction for a diagonal profile in the Onilus Region:
	lon = 20 + np.arange(11)
	lat = 40 + np.arange(11)
	#d = np.sqrt( x**2 + y**2)
	z = extract_MOLA_profile(lat, lon)
	print(z)
