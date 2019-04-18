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
	"""Class to read, store, and process seisware horizons,
		specifically for the SWIM project"""

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

	def estimate_epsilon_along_track_MOLA_minima(self, sub_horiz, surf_horiz, MOLA_file, filepath='./', orbit_list=[]):
		"""Function to estimate the real dielectric constant 
		for subsurface reflectors using the MOLA minima
		extrapolation method developed in MATLAB by Ali 
		Bramson. Uses a function called from general_functions.py
		to do the majority of the heavy lifting.

		:param str sub_horiz: subsurface horizon/reflector to
			estimate the real dielectric constant for
		:param str surf_horiz: surface horizon/reflector
		:param str MOLA_file: filepath to MOLA geotiff
			file, used to extract MOLA elevation
		:param str filepath: filepath to export outputs
		"""

		# Check if filepaths exist, and if not make them
		figpath = filepath + 'Plots/'
		if not os.path.exists(filepath):
			os.makedirs(filepath)
		if not os.path.exists(figpath):
			os.makedirs(figpath)
		
		# Start writing file for the segments summary:
		seg_sum_outpath = filepath+sub_horiz + '_summary.csv'
		seg_sum_header = ['Orbit','Center Lat','Center Lon','Median Eps','Mean Eps','Std Eps','IQR Eps','Min1 Lon','Min1 Lat','Min2 Lon','Min2 Lat']
		fout = open(seg_sum_outpath, 'wb')
		seg_sum_fout = csv.writer(fout)
		seg_sum_fout.writerow(seg_sum_header)
		
		# Initialize 

		# Find the orbits which contain the subsurface
		#	horizon:
		TWT_sub = self.horizons[sub_horiz]
		TWT_surf = self.horizons[surf_horiz]
		isind_sub = np.where( (TWT_sub != -9999.99) & (TWT_sub != np.nan) )
		if orbit_list == []:
			orbit_list = np.unique( self.orbit[isind_sub] )
		# Cycle through those orbits and do the dielectric
		# 	estimation:
		OO = np.size(orbit_list)
		for oo in range(OO):
			orbind = np.where( (self.orbit == orbit_list[oo]) & (TWT_surf != -9999.99) )
			epsilon, z_sub, seg_sum = estimate_epsilon_MOLA_minima_extrapolation(orbit_list[oo], self.trace[orbind], self.lat[orbind], self.lon[orbind], TWT_surf[orbind], TWT_sub[orbind], MOLA_file, figpath)
			if np.size(seg_sum)==11:
				seg_sum = seg_sum.flatten()
				seg_sum_fout.writerow(seg_sum)
			if np.size(seg_sum)>11:
				seg_sum_fout.writerows(seg_sum)
			if oo == 0:
				orball = self.orbit[orbind]
				traceall = self.trace[orbind]
				lonall = self.lon[orbind]
				latall = self.lat[orbind]
				epsilonall = epsilon
				z_suball = z_sub
				TWT_surfall = TWT_surf[orbind]
				TWT_suball = TWT_sub[orbind]
			else:
				orball = np.concatenate( (orball, self.orbit[orbind]), axis=None)
				traceall = np.concatenate( (traceall, self.trace[orbind]), axis=None)
				lonall = np.concatenate( (lonall, self.lon[orbind]), axis=None)
				latall = np.concatenate( (latall, self.lat[orbind]), axis=None)
				epsilonall = np.concatenate( (epsilonall, epsilon), axis=None)
				z_suball = np.concatenate( (z_suball, z_sub), axis=None)
				TWT_surfall = np.concatenate( (TWT_surfall, TWT_surf[orbind]), axis=None)
				TWT_suball = np.concatenate( (TWT_suball, TWT_sub[orbind]), axis=None)
		# Median epsilon of whole dataset:
		med_eps_all = np.nanmedian(epsilonall)
		print('Median Epsilon for all reflectors = {}').format(med_eps_all)
		# Depth correct using median epsilon:
		depth_med_eps = depth_correct_radar(TWT_surfall, TWT_suball, med_eps_all)

		# Prepare dataset for output:
		ii = np.where( (TWT_suball != -9999.99) & (TWT_suball != np.nan))
		outdat = np.column_stack((orball[ii], traceall[ii], latall[ii], lonall[ii], epsilonall[ii], depth_med_eps[ii], z_suball[ii], TWT_surfall[ii], TWT_suball[ii]))
		outhead = ['Orbit', 'Trace', 'Latitude', 'Longitude', 'Epsilon', 'Depth (median eps ='+str(med_eps_all)+')','Z_Sub','TWT surf (ns)','TWT sub (ns)']
		write_csv_file(filepath+sub_horiz+'_result.csv', outdat, outhead)

if __name__ == '__main__':

	MOLAfile='/Users/eric/Documents/orig/supl/MOLA/DEM_global_mola128ppd_merged_mola64ppd/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif'
	
	datafile = '../Horizon_Export/2019_04_08.txt'
	filepath = '../Horizon_Export/'

	data = swim_horizons(datafile, onilus)
#	data.write_horizon_csvs(filepath)
	#data.depth_correct('lda_sub1_EP','lda_surf_EP',3, filepath)
	#data.estimate_epsilon_targ_horiz('plains_sub1_EP','surf_EP','target_base_EP', filepath)

	# Test MOLA minima method:
	orbind = np.where((data.orbit == 3550401) & (data.horizons['surf_EP'] != -9999.99))
	epsilon, z_sub, seg_sum = estimate_epsilon_MOLA_minima_extrapolation(3550401, data.trace[orbind], data.lat[orbind], data.lon[orbind], data.horizons['surf_EP'][orbind], data.horizons['plains_sub1_EP'][orbind], MOLAfile, '../Ali_Dielectric_Code/Plots_SHARAD_Estimations/')
	print seg_sum
