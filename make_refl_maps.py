import numpy as np
import gdal
import osgeo.ogr, osgeo.osr
import matplotlib.pyplot as plt
import os
import sys 
import csv
from mars_projections import *
from seisware_horizons import *
#import gdal_merge
#import subprocess

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

class swim_horizons(seisware_horizons):
	"""Class to read, store, and process seisware horizons"""

	def __init__(self, datafile, src, dst=mars_2000):
		"""Initiate based on seisware horizons class

		:param datafile: seisware output to be read in
		:param src: Seisware project projection
		:param dst: Desired output projection, by default
			Mars 2000
		"""

		super(swim_horizons, self).__init__(datafile, src, dst)

	def make_conf_refl(self, sub_horiz, res, outfile):
		"""Produces confidence map based on the presence
		or non-presence of subsurface radar reflectors

		:param str sub_horiz: list of names of mapped 
			subsurface horizons; must match options
			found in self.horiz_names
		:param str outfile: file to write results to
		"""

		print 
		print('Producing tif file of reflection confidence value from the following horizons')
		for horiz in sub_horiz:
			print('    {}'.format(horiz))
		print

		conf_refl = np.zeros(self.num_points)
		num_sub_horiz = len(sub_horiz)
		for N in range(num_sub_horiz):
			temp_horiz = self.horizons[sub_horiz[N]]
			mapped_refl = np.argwhere(temp_horiz != self.nan_const)
			conf_refl[mapped_refl] = 1
		
		self.conf_refl = conf_refl

		convert_xyz_to_raster(self.lon, self.lat, conf_refl, res, 'maximum', outfile)

if __name__ == '__main__':
	# Usage: python SWIM_process_picks pick_file region sub_horiz
	# Inputs:
	res = 0.1 #Resolution in degrees lon/lat
	if len(sys.argv) < 3:
		print
		print 'ERROR: Incorrect number of arguments (2 or more expected)'
		print
		print 'Usage: "python SWIM_process_picks pick_file region sub_horiz1, sub_horiz2 ...'
		print 'pick_file = seisware horizons output file'
		print 'region = onilus, acidalia, utopia, or arcadia'
		print '(Optional): sub_horiz1, 2, etc = names of subsurface horizons for '
		print '			producing Reflection Confidence Map'
		print
		sys.exit(1)

	pick_file = sys.argv[1]
	region = sys.argv[2]
	num_sub_horiz = len(sys.argv) - 3
	sub_horiz = []
	if num_sub_horiz > 0:
		for N in range(num_sub_horiz):
			sub_horiz.append(sys.argv[3+N])

	if region == 'onilus': region = onilus
	if region == 'acidalia': region = acidalia
	if region == 'utopia': region = utopia
	if region == 'arcadia': region = arcadia

	pick_path = os.path.dirname(pick_file) + '/' 

	horizons = swim_horizons(pick_file, region) # Read in Horizons
	horizons.write_horizon_csvs(pick_path) # Save csv files of individual horizons
	if sub_horiz == 'all':
		sub_horiz = horizons.horiz_names
	if sub_horiz:
		outfile = pick_file[:-4] + '_refl_conf.tif'
		horizons.make_conf_refl(sub_horiz, res, outfile)
