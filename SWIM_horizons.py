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
from general_functions import *


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


	def estimate_epsilon_targ_horiz(self, sub_horiz, surf_horiz, targ_horiz, filepath='./'):
		"""Estimates epsilon for a subsurface horizon by
		prescribing a target horizon, then uses the median dielectric
		constant found by this method to depth correct the sub_horiz

		:param str sub_horiz: subsurface horizon to depth-correct
		:param str surf_horiz: surface horizon
		:param str targ_horiz: target horizon to depth-correct to
		:param str filepath: option filepath where depth data is saved 
		"""

		sub = self.horizons[sub_horiz]
		surf = self.horizons[surf_horiz]
		targ = self.horizons[targ_horiz]

		ind = (sub != self.nan_const) & (surf != self.nan_const) & (targ != self.nan_const)

		epsilon = estimate_epsilon_targetTWT(surf[ind], sub[ind], targ[ind])
		conf = convert_epsilon_to_conf(epsilon)
		eps_med = np.median(epsilon)
		eps_std = np.std(epsilon)
		print
		print('Estimating epsilon for horizon {}'.format(sub_horiz))
		print('Median Epsilon = {}'.format(eps_med))
		print('StDev = {}'.format(eps_std))

		depth = depth_correct_radar(surf[ind], sub[ind], eps_med)
		
		outfile = filepath + sub_horiz + '_eps.csv'

		print('Saving estimated epsilon and depths to {}'.format(outfile))	
		print('Number of Radargrams = {}'.format(len(set(self.orbit[ind]))))
		print('Number of Points = {}'.format(len(self.orbit[ind])))
		print('Median Depth = {}'.format(np.median(depth)))


		out_header = ['Orbit', 'Trace', 'Latitude', 'Longitude', 'Epsilon', 'Cons Value', 'Depth (eps='+str(eps_med)+')', 'TWT surf (ns)', 'TWT sub (ns)', 'TWT targ (ns)']
		out_data = np.column_stack((self.orbit[ind], self.trace[ind], self.lat[ind], self.lon[ind], epsilon, conf, depth, surf[ind], sub[ind], targ[ind]))

		write_csv_file(outfile, out_data, out_header)


if __name__ == '__main__':
	
	datafile = '../Horizon_Export/2019_02_22.txt'
	filepath = '../Horizon_Export/'

	data = swim_horizons(datafile, onilus)
#	data.write_horizon_csvs(filepath)
	#data.depth_correct('lda_sub1_EP','lda_surf_EP',3, filepath)
	data.estimate_epsilon_targ_horiz('plains_sub1_EP','surf_EP','target_base_EP', filepath)
