import numpy as np
import osgeo.ogr, osgeo.osr
import matplotlib.pyplot as plt
import os 
import csv
from mars_projections import *
from general_functions import *


class seisware_horizons(object):
	"""This class holds information on seisware horizons"""
	
	def __init__(self, datafile, src, dst=mars_2000):
		"""Read in tab-delimited seisware horizons
		export and return the data contained within.
		Also converts the coordinates of the horizons,
		by default from the project Albers projection
		to the standard Mars 2000

		:param datafile: Datafile to be read in
		:param src: Seisware Project Projection
		:param dst: Desired output projection, by 
			default Mars 2000
		"""

		# Read in data file:
		orbit = []
		trace = []
		X = []
		Y = []
		horiz_name = []
		picks = {} # Pick information will be stored in a dictionary
		f = open(datafile, 'r')
		with f:
			content = f.readlines()
			line_number = 0
			for line in content:
				row = line.split()
				if line_number == 0:
					n_horiz = int(row[0]) # Number of Horizons
					nan_const = np.array(row[1], dtype='float64') # NaN Constant
				if (line_number>0 and line_number<=n_horiz):
					horiz_name.append(row[0])
					picks[row[0]] = []
				if line_number > n_horiz: 
					orbit.append(int(row[0]))
					trace.append(int(row[1]))
					X.append(row[2])
					Y.append(row[3])
					for i in range(n_horiz):
						picks[horiz_name[i]].append(row[4+i])
				line_number = line_number+1
	
		# Convert X, Y Projected Points to Lat, Lon in
		#	Mars 2000 Reference System:
		X = np.array(X, dtype = 'float')
		Y = np.array(Y, dtype = 'float')
		points = []
		for xx in range(len(X)):
			points.append([X[xx], Y[xx]])
	
		lonlat = transform_points(points, src, dst)	
		lon = []
		lat = []
		for xx in range(len(lonlat)):
			point = lonlat[xx]
			lat.append(point[1])
			lon.append(point[0])
		lat = np.array(lat)
		lon = np.array(lon)

		# Convert horizon picks from "seisware units" to 
		# 	nanoseconds:
		seis_conversion = 100 # 100 ns/seisware unit
		for i in range(n_horiz):
			picks[horiz_name[i]] = np.array(picks[horiz_name[i]], dtype='float64') # Convert to float
			ind = np.argwhere(picks[horiz_name[i]] != nan_const) # Non-nan index
			picks[horiz_name[i]][ind] = picks[horiz_name[i]][ind] * seis_conversion # convert to ns where non-nan

		# Populate object:
		self.num_horizons = n_horiz
		self.nan_const = nan_const
		self.horiz_names = horiz_name
		self.lon = lon
		self.lat = lat
		self.orbit = np.array(orbit)
		self.trace = np.array(trace)
		self.X = X
		self.Y = Y
		self.horizons = picks
		self.num_points = len(self.X)
	
		print()
		print('Read in Seisware File {}'.format(datafile))
		print('Number of Radargrams: {}'.format(len(set(self.orbit))))
		print('Number of Points: {}'.format(self.num_points))
		print('{} Horizons include:'.format(self.num_horizons))
		for horiz in self.horiz_names:
			print('    {}'.format(horiz))
		print('Longitude from {0} to {1}'.format(min(lon), max(lon)))
		print('Latitude from {0} to {1}'.format(min(lat), max(lat)))

	def write_horizon_csvs(self, filepath='./'):
		"""Reads in seisware class recording horizons
		and outputs csv file for each containing the info,
		with the coordinates converted to Mars 2000
		Reference system

		:param str filepath: optional filepath indicating
			where to save the horizon files; default is
			'./' for current folder.
		"""
		print()
		print('Saving Horizon Files')
		for horiz in self.horiz_names:
			# Index to select only non-nan data:
			ind = (self.horizons[horiz] != self.nan_const)
			out_data = np.column_stack((self.orbit[ind], self.trace[ind], self.lat[ind], self.lon[ind], self.horizons[horiz][ind]))
			out_header = ['Orbit', 'Trace', 'Latitude', 'Longitude', 'TWT (ns)']
			outfile = filepath + horiz + '.csv' # Filename
			# Export to csv file:
			print('    Saving horizon {0} as {1}'.format(horiz, outfile))
			print('        Contains {} points.'.format(len(self.orbit[ind])))
			write_csv_file(outfile, out_data, out_header)
	def depth_correct(self, sub_horiz, surf_horiz, epsilon, filepath='./'):
		"""Depth-correct given subsurface horizon
		using the specified surface horizon and a specified
		dielectric constant value epsilon

		:param str sub_horiz: subsurface horizon to depth-correct
		:param str surf_horiz: surface horizon
		:param float epsilon: dielectric constant value
		:param str filepath: optional filepath where depth data is 
			saved
		"""

		sub = self.horizons[sub_horiz]
		surf = self.horizons[surf_horiz]

		ind = ( sub != self.nan_const ) & ( surf != self.nan_const )
		
		depth = depth_correct_radar(surf[ind], sub[ind], epsilon)

		outfile = filepath + sub_horiz + '_depth_' + str(epsilon) + '.csv'

		print()
		print('Depth-corrected {0} using epsilon = {1}'.format(sub_horiz, epsilon))
		print('Median depth = {}'.format(np.median(depth)))
		print('Saving as {}'.format(outfile))
		print('Number of Radargrams = {}'.format(len(set(self.orbit[ind]))))
		print('Number of Points = {}'.format(len(self.orbit[ind])))
		print()

		out_header = ['Orbit', 'Trace', 'Latitude', 'Longitude', 'Depth', 'TWT surf (ns)', 'TWT sub (ns)']
		out_data = np.column_stack((self.orbit[ind], self.trace[ind], self.lat[ind], self.lon[ind], depth, surf[ind], sub[ind]))

		write_csv_file(outfile, out_data, out_header)


