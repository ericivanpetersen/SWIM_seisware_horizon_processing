
import numpy as np
import gdal
import csv

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

