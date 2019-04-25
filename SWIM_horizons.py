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
		print('Number of Orbits = {}').format(OO)
		# "Bookmark" setting: start from specific orbit
		orb_success = 0
		while orb_success == 0:
			start_orb = input('Orbit start number? (if start from beginning, enter 0): ')
			if start_orb == 0:
				oo1 = 0
				seg_sum_outpath = filepath+sub_horiz + '_summary.csv'
				result_outpath = filepath+sub_horiz + '_results.csv'
				depth_result_outpath = filepath+sub_horiz + '_depth_results.csv'
				orb_success = 1
			else:
				oo1 = np.argwhere( orbit_list == start_orb )
				if np.size(oo1)==0:
					print("Orbit does not exist; try again")
				else:
					seg_sum_outpath = filepath+sub_horiz + '_summary'+str(oo1.flatten())+'.csv'
					result_outpath = filepath+sub_horiz + '_results'+str(oo1.flatten())+'.csv'
					depth_result_outpath = filepath+sub_horiz + '_depth_results'+str(oo1.flatten())+'.csv'
					orb_success = 1
		OO = np.size(orbit_list)
		# Start writing file for the segments summary:
		seg_sum_header = ['Orbit','Center Lat','Center Lon','Median Eps','Mean Eps','Std Eps','IQR Eps','Min1 Lat','Min1 Lon','Min2 Lat','Min2 Lon']
		result_header = ['Orbit','Trace','Latitude','Longitude','Epsilon','Cons','Z_Sub','TWT_surf (ns)','TWT_sub (ns)']
		fout1 = open(seg_sum_outpath, 'wb')
		fout2 = open(result_outpath,'wb')
		seg_sum_fout = csv.writer(fout1)
		seg_sum_fout.writerow(seg_sum_header)
		result_fout = csv.writer(fout2)
		result_fout.writerow(result_header)
		for oo in range(oo1,OO):
			orbind = np.where( (self.orbit == orbit_list[oo]) & (TWT_surf != -9999.99) )
			epsilon, z_sub, seg_sum = estimate_epsilon_MOLA_minima_extrapolation(orbit_list[oo], self.trace[orbind], self.lat[orbind], self.lon[orbind], TWT_surf[orbind], TWT_sub[orbind], MOLA_file, figpath)
			# Write data to running results file:
			if np.size(seg_sum) > 0:
				lind = np.invert(np.isnan(epsilon))
				cons = convert_epsilon_to_conf(epsilon[lind])
				result_rows_out = np.column_stack([self.orbit[orbind][lind], self.trace[orbind][lind], self.lat[orbind][lind], self.lon[orbind][lind], epsilon[lind], cons, z_sub[lind], TWT_surf[orbind][lind], TWT_sub[orbind][lind]])
				result_fout.writerows(result_rows_out) 
			# Write data to running summary file:
			if np.size(seg_sum)==11:
				seg_sum = seg_sum.flatten()
				seg_sum_fout.writerow(seg_sum)
			if np.size(seg_sum)>11:
				seg_sum_fout.writerows(seg_sum)
			# Concatenate data for final depth results output:
			if oo == oo1:
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
		print('Std Epsilon = {}').format(np.nanstd(epsilonall))
		# Depth correct using median epsilon:
		depth_med_eps = depth_correct_radar(TWT_surfall, TWT_suball, med_eps_all)

		# Prepare dataset for output:
		ii = np.where( (TWT_suball != -9999.99) & (TWT_suball != np.nan))
		outdat = np.column_stack((orball[ii], traceall[ii], latall[ii], lonall[ii], epsilonall[ii], depth_med_eps[ii], z_suball[ii], TWT_surfall[ii], TWT_suball[ii]))
		outhead = ['Orbit', 'Trace', 'Latitude', 'Longitude', 'Epsilon', 'Depth (median eps ='+str(med_eps_all)+')','Z_Sub','TWT surf (ns)','TWT sub (ns)']
		write_csv_file(filepath+sub_horiz+'_depth_result.csv', outdat, outhead)

	def estimate_epsilon_nearest_plains_elevation(self, sub_horiz, surf_horiz, plains_horiz, filepath):
		"""Function to estimate real dielectric constant,
		using a method that  Eric Petersen developed for fast, simple
		dielectric estimation for LDAs in Deuteronilus Mensae.
		This function is a relatively faithful reproduction of 
		the MATLAB script used to estimate dielectric constant 
		for the GRL paper Petersen et al, 2018.
		Script finds the nearest plains elevation based on
		the picked horizons and corrects
		the subsurface reflector to align with it. Surface 
		points are not calibrated to MOLA, to avoid issues
		with off-nadir LDA returns not being well represented
		by nadir MOLA.

		:param str sub_horiz: subsurface horizon
		:param str surf_horiz: surface horizon over subsurface horizon
		:param str plains_horiz: plains horizon located adjacent
			to LDA/other feature with subsurface horizon
		"""

		if not os.path.exists(filepath):
			os.makedirs(filepath)
		prox_lim = 1500
		# Prepare to save results:
		result_outpath = filepath + sub_horiz + '_results.csv'
		seg_sum_outpath = filepath + sub_horiz + '_summary.csv'
		result_header = ['Orbit','Trace','Latitude','Longitude','Epsilon','Cons','Z_Sub','TWT_surf (ns)','TWT_sub (ns)']
		seg_sum_header = ['Orbit','Center Lat','Center Lon','Median Eps','Mean Eps','Std Eps','IQR Eps','Lat Plains','Lon Plains','Z_Plains']
		fout1 = open(seg_sum_outpath, 'wb')
		fout2 = open(result_outpath, 'wb')
		seg_sum_fout = csv.writer(fout1)
		result_fout = csv.writer(fout2)
		seg_sum_fout.writerow(seg_sum_header)
		result_fout.writerow(result_header)
		# Find the orbits with sub_horizon and plains_horizon present:
		TWT_sub = self.horizons[sub_horiz]
		TWT_surf = self.horizons[surf_horiz]
		TWT_plains = self.horizons[plains_horiz]
		orb_sub = self.orbit[np.where( (TWT_sub != self.nan_const))]
		orb_plains = self.orbit[np.where( (TWT_plains != self.nan_const))]
		orb_surf = self.orbit[np.where( (TWT_surf != self.nan_const))]
		orbit_list = np.intersect1d(orb_sub, orb_plains)
		orbit_list = np.intersect1d(orbit_list, orb_surf)
		OO = np.size(orbit_list)
		print('Number of Orbits = {}').format(OO)
		# Loop through orbits:
		oo = 0
		epsall = []
		while oo < OO:
			orb_ind = np.where( self.orbit == orbit_list[oo])
			z_surf = -0.3 * TWT_surf[orb_ind]/2
			z_surf_plains = -0.3 * TWT_plains[orb_ind]/2 
			z_surf = z_surf.flatten()
			sub_ind = np.where( TWT_sub[orb_ind] != self.nan_const)
			sub_ind = sub_ind[0]
			tdiff = np.diff(self.trace[orb_ind][sub_ind]) # trace-length between mapped subsurface horizons
			TT = np.size(sub_ind) # number of traces with subsurface horizons
			tt = 0
			while tt < TT:
				# Find unbroken segment of reflector:
				tindend = np.where(tdiff[tt:] != 1)
				tindend = tindend[0]
				if np.size(tindend) == 0:
					trace_ind = np.arange(tt,TT)
					tt = TT
				else:
					trace_ind = np.arange(tt,tindend[0]+tt)
					tt = tindend[0]+tt+1
				trace_ind = sub_ind[trace_ind]
				if np.size(trace_ind) > 0:
					# Find nearest plains surface which is below elevation of LDA:
					z_lda_srf_min = np.amin( z_surf[trace_ind] ) # Minimum elevation of LDA surf
					p_ind1 = np.where( (z_surf_plains < z_lda_srf_min) & ( TWT_plains[orb_ind] != self.nan_const) & (self.trace[orb_ind] < self.trace[orb_ind][trace_ind][0])) #Where the plains reflector mapped at lower elevation than the LDA surface to the left 
					p_ind2 = np.where( (z_surf_plains < z_lda_srf_min) & (TWT_plains[orb_ind] != self.nan_const) & (self.trace[orb_ind] > self.trace[orb_ind][trace_ind][-1])) #Where the plains reflector mapped at lower elevation than the LDA surface to the right
					p_ind1 = p_ind1[0]
					p_ind2 = p_ind2[0]
					# Special Case for reflectors found near upper plains unit: default to plains near the south for dielectric estimation				
					if (np.mean( self.lon[orb_ind][trace_ind]) > 25.5) & (np.mean( self.lon[orb_ind][trace_ind]) < 30.5) & (np.mean( self.lat[orb_ind][trace_ind]) > 43.5):
						p_ind2 = np.where( (z_surf_plains < z_lda_srf_min) & (TWT_plains[orb_ind] != self.nan_const) & (self.lat[orb_ind] < 43.5))
						p_ind2 = p_ind2[0]
						p_ind1 = []
					
					if (np.size(p_ind1) == 0) & (np.size(p_ind2) == 0):
						continue
					elif (np.size(p_ind2) == 0):
						side = 1 # to the left, to the left
						prox1 = self.trace[orb_ind][trace_ind][0] - self.trace[orb_ind][p_ind1][-1]
					elif (np.size(p_ind1) == 0):
						side = 2 # to the right, to the right
						prox2 = self.trace[orb_ind][p_ind2][0] - self.trace[orb_ind][trace_ind][-1]
					else:
						prox1 = self.trace[orb_ind][trace_ind][0] - self.trace[orb_ind][p_ind1][-1]
						prox2 = self.trace[orb_ind][p_ind2][0] - self.trace[orb_ind][trace_ind][-1]
						if prox1 < prox2:
							side = 1
						else:
							side = 2
					if side == 1:
						prox = prox1 
						z_p = z_surf_plains[p_ind1[-1]]
						lon_p = self.lon[orb_ind][p_ind1[-1]]
						lat_p = self.lat[orb_ind][p_ind1[-1]]
					else:
						prox = prox2
						z_p = z_surf_plains[p_ind2[0]]
						lon_p = self.lon[orb_ind][p_ind2[0]]
						lat_p = self.lat[orb_ind][p_ind2[0]]
					if prox > prox_lim:
						continue					

					depth = z_surf[trace_ind] - z_p
					depth = depth.flatten()
					eps = estimate_epsilon( TWT_surf[orb_ind][trace_ind], TWT_sub[orb_ind][trace_ind], depth)
					cons = convert_epsilon_to_conf(eps)
					# Save results summary:
					if (np.nanmedian(eps) < 12) & (np.nanmedian(eps)>1) & (np.nanstd(eps)<10):
						print('Line {0}; Epsilon = {1}').format(orbit_list[oo],np.nanmedian(eps))
						epsall = np.concatenate( (epsall, eps), axis=None)
						center_lon = np.mean(self.lon[orb_ind][trace_ind])
						center_lat = np.mean(self.lat[orb_ind][trace_ind])
						seg_sum = [orbit_list[oo], center_lat, center_lon, np.nanmedian(eps), np.nanmean(eps), np.nanstd(eps), np.subtract(*np.nanpercentile(eps, [75, 25])), lat_p, lon_p, z_p]
						seg_sum_fout.writerow(seg_sum)

						# Save results:	
						seg_res = np.column_stack([self.orbit[orb_ind][trace_ind], self.trace[orb_ind][trace_ind], self.lat[orb_ind][trace_ind], self.lon[orb_ind][trace_ind], eps, cons, z_surf[trace_ind]-depth, TWT_surf[orb_ind][trace_ind], TWT_sub[orb_ind][trace_ind]])
						result_fout.writerows(seg_res)
			oo = oo+1
		print('Median Epsilon = {}').format(np.median(epsall)) 
		print('Mean Epsilon = {}').format(np.mean(epsall))
		print('Stdev Epsilon = {}').format(np.std(epsall))
		print('Max Epsilon = {}').format(np.amax(epsall))
				
if __name__ == '__main__':

	MOLAfile='/Users/eric/Documents/orig/supl/MOLA/DEM_global_mola128ppd_merged_mola64ppd/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif'
	
	datafile = '../Horizon_Export/2019_04_08.txt'
	filepath = '../Horizon_Export/Dielectric_NPE/'

	data = swim_horizons(datafile, onilus)
	data.estimate_epsilon_nearest_plains_elevation('lda_sub1_EP','lda_surf_EP','surf_EP',filepath)
	data.depth_correct('lda_sub1_EP','lda_surf_EP',3, filepath)
	data.depth_correct('up_sub1_EP','surf_EP',4,filepath)
	data.depth_correct('plains_sub1_EP','surf_EP',3.7,filepath)
	data.depth_correct('plains_sub1_EP','surf_EP',5.2,filepath)
