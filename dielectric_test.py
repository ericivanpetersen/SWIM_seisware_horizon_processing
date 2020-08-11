import numpy as np
import pandas as pd
from general_functions import *
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

MOLA_file='/Users/eric/Documents/orig/supl/MOLA/DEM_global_mola128ppd_merged_mola64ppd/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif'

def generate_random_profile(LL, h, z0):
	"""
	:param int L: length of profile
	:param float h: width of normal distribution
		for producing random profile
	:param float z0: seed value for initial z
	:output float Z: Z value for randomly
		generated profile
	"""

	Z = np.full(LL, 0.0)
	for ll in range(LL):
		if ll == 0:
			Z[ll] = z0 + h * np.random.randn()
		if ll > 0:
			Z[ll] = Z[ll-1] + h * np.random.randn()

	return Z

def monte_carlo_eps_estimate(delta_TWT, Z_surf, z_base_0, h_rms, lead_in, num_iter=10000):
	"""
	Experiment to use monte carlo generation of 
	basal reflectors to estimate the dielectric 
	constant as well as basal topography of the
	subsurface reflector
	
	:param float delta_TWT: difference in TWT 
		between surface & subsurface
	:param float Z_surf: elevation of the surface
	:param float z_base_0: seed value for creating 
		synthetic profiles
	:param int num_iter: number of monte carlo
		iterations.
	"""
	
	prof_length = np.size(delta_TWT) # length of the subsurface profile

	# arrays to fill with data from the experiment:
	eps_est = np.full(num_iter, 0.0) # best-fit epsilon for each profile
	base_MAE = np.full(num_iter, 0.0) # MAE between data and synthetic profile
	base_R2 = np.full(num_iter, 0.0) # R2 of fit between data and synthetic profile
	base_prof = np.full([num_iter, prof_length], z_base_0) # generated profiles

	for nn in range(num_iter):
		if lead_in > 0:
			z_0_nn = z_base_0
			for pp in range(lead_in):
				z_0_nn = z_0_nn + h_rms * np.random.randn()
		temp_base = generate_random_profile(prof_length, h_rms, z_0_nn)
		base_prof[nn, :] = temp_base
		H = Z_surf - temp_base
		H = H.reshape((-1, 1))
		linmod = LinearRegression(fit_intercept=False).fit(H, delta_TWT)
		eps_est[nn] = (linmod.coef_[0] * 0.3/2) ** 2
		base_R2[nn] = linmod.score(H, delta_TWT)
		base_MAE[nn] = mean_absolute_error(H, (0.3 / np.sqrt(eps_est[nn]) * delta_TWT/2))

	max_R2 = np.amax(base_R2)
	mean_eps = np.mean(eps_est)
	std_eps = np.std(eps_est)
	min_MAE = np.amin(base_MAE[np.argwhere( (eps_est > mean_eps-std_eps) & (eps_est < mean_eps+std_eps))])
	#best_ind = np.argwhere(base_R2 == max_R2)[0]
	#print(best_ind)
	#print("Best fit profile eps = {0} at R2 = {1}".format(eps_est[best_ind], max_R2))
	best_ind = np.argwhere(base_MAE == min_MAE)[0]	
	print("Best fit profile eps = {0} at MAE = {1}".format(eps_est[best_ind], min_MAE))
	print("Median Epsilon = {}".format(np.median(eps_est)))
	print("Mean Epsilon = {0} +- {1}".format(np.mean(eps_est), np.std(eps_est)/np.sqrt(np.size(eps_est))))
	print("Std Eps = {}".format(np.std(eps_est)))
	best_prof = base_prof[best_ind, :]
 
	return eps_est, base_MAE, base_prof, best_ind

def main():
	sub_dat = pd.read_csv('../../Horizon_Export/lda_sub1_EP_depth_3.csv')
	orbit = sub_dat['Orbit'].to_numpy()
	trace = sub_dat['Trace'].to_numpy()
	orbind = np.where((orbit == 1736801) & (trace>335) & (trace<450))[0]
#	orbind = np.where((orbit == 755102) & (trace > 550) & (trace< 700))
#	orbind = np.where((orbit == 753802) & (trace > 700) & (trace<800))[0]
#	orbind = np.where((orbit == 690502) & (trace > 900) & (trace < 1000))
#	orbind = np.where((orbit == 1274601) & (trace > 1140) & (trace < 1300))
	sub_TWT = sub_dat['TWT sub (ns)'].to_numpy()[orbind]
	surf_TWT = sub_dat['TWT surf (ns)'].to_numpy()[orbind]
	trace = sub_dat['Trace'].to_numpy()[orbind]
	orbit = sub_dat['Orbit'].to_numpy()[orbind]
	lat = sub_dat['Latitude'].to_numpy()[orbind]
	lon = sub_dat['Longitude'].to_numpy()[orbind]
	z_mola = extract_MOLA_profile(lon, lat, MOLA_file)
	z_mola = np.squeeze(z_mola)

	plains_dat = pd.read_csv('../../Horizon_Export/surf_EP.csv')
	orb_pl = plains_dat['Orbit']
	trace_pl = plains_dat['Trace']
#	pind = np.where((orb_pl == 1736801) & (trace_pl > 349) & (trace_pl < 372) )
	pind = np.where((orb_pl == 1736801) & (trace_pl > 335) & (trace_pl < 400) )
#	pind = np.where((orb_pl == 755102) & (trace_pl > 520) & (trace_pl < 600) )
#	pind = np.where((orb_pl == 753802) & (trace_pl > 650) & (trace_pl < 750))
#	pind = np.where((orb_pl == 690502) & (trace_pl > 850) & (trace_pl < 950))
#	pind = np.where((orb_pl == 1274601) & (trace_pl > 1100) & (trace_pl < 1200))
	plains_TWT = plains_dat['TWT (ns)'].to_numpy()[pind]
	plains_trace = plains_dat['Trace'].to_numpy()[pind]
	plains_lon = plains_dat['Longitude'].to_numpy()[pind]
	plains_lat = plains_dat['Latitude'].to_numpy()[pind]
	plains_z_mola = extract_MOLA_profile(plains_lon, plains_lat, MOLA_file)
	plains_z_mola = np.squeeze(plains_z_mola)
	#[eps_best, depth_best, r2] = estimate_epsilon_min_R2(trace, surf_TWT, sub_TWT, plot=True) 

	delta_TWT = sub_TWT - surf_TWT 
	Z_surf = - surf_TWT * 0.3/2
	plains_z = - plains_TWT * 0.3/2 # Elevation below MRO

	MOLA_conversion = np.mean(plains_z_mola) - np.mean(plains_z)
	plains_z = plains_z + MOLA_conversion
	Z_surf = Z_surf + MOLA_conversion

#	Z_surf = z_mola
#	plains_z = plains_z_mola
	plains_z_diff = plains_z[1:] - plains_z[:-1]
	h_rms = np.sqrt(np.mean(plains_z_diff ** 2))
	z_0 = np.mean(plains_z[-4:-1])
	z_0_mola = np.mean(plains_z_mola[-4:-1])
	m = np.mean(plains_z_diff)
	lead_in = np.abs(trace[0] - plains_trace[-1]) 
	print("H_rms = {}".format(h_rms))
	print("z0 = {}".format(z_0))
	print("m = {}".format(m))
	print("Lead-in = {}".format(lead_in))
	
	# Nearest Plains; MOLA vs Raw SHARAD:
	eps_best_so = estimate_epsilon(surf_TWT, sub_TWT, (Z_surf - z_0))
	print("Mean Eps SHARAD-Only: {0}  STD: {1}".format(np.mean(eps_best_so), np.std(eps_best_so)))
	eps_best_MOLA = estimate_epsilon(surf_TWT, sub_TWT, (z_mola - z_0_mola))
	print("Mean Eps, MOLA-fixed: {0}  STD: {1}".format(np.mean(eps_best_MOLA), np.std(eps_best_MOLA)))	

	# Plot comparison between MOLA & TWT raw:
	plt.plot(trace, Z_surf, 'k-')
	plt.plot(plains_trace, plains_z,'k-')
	plt.plot(trace, z_mola,'r-')
	plt.plot(plains_trace, plains_z_mola,'r-')
	plt.plot(trace, Z_surf - delta_TWT * 0.3/2, 'k-')
	plt.plot(trace, np.full(np.size(trace), z_0), 'b--')
	plt.xlabel('Trace')
	plt.ylabel('MOLA Elevation')
	plt.show()
	plt.plot(trace, eps_best_so, 'k.')
	plt.plot(trace, eps_best_MOLA, 'ro')
	plt.legend(['SHARAD-Only', 'MOLA'])
	plt.xlabel('Trace')
	plt.ylabel('Epsilon')
	plt.show()

	plt.subplot(223)
	plt.plot(-surf_TWT, z_mola, 'k.')
	plt.plot(-plains_TWT, plains_z_mola, 'ko')
	plt.legend('LDA Surface', 'Plains Surface')
	plt.xlabel('SHARAD-Only')
	plt.ylabel('MOLA')
	plt.subplot(224)
	plt.plot(eps_best_so, eps_best_MOLA, 'k.')
	plt.xlabel('SHARAD-Only')
	plt.ylabel('MOLA')
	plt.show()
	
	eps_est, base_MAE, profs, best_ind = monte_carlo_eps_estimate(delta_TWT, Z_surf, z_0, h_rms, lead_in)

	print(np.size(profs))

	plt.hist(eps_est, 50)
	plt.xlabel('Epsilon')
	plt.show()

	plt.hist(base_MAE, 50)
	plt.xlabel('MAE')
	plt.show()

	for n in range(5):
		plt.plot(trace, Z_surf)
		plt.plot(trace, profs[n,:].flatten())
		plt.plot(trace, Z_surf - 0.3 /np.sqrt(eps_est[n]) * delta_TWT/2, 'k--')
		plt.plot(trace, Z_surf - delta_TWT * 0.3/2)
		plt.plot(plains_trace, plains_z)
		plt.legend(['Surf','Generated Base','Depth-Corrected Reflector','Uncorrected Reflector','Plains'])
		plt.xlabel('Trace')
		plt.ylabel('meters below MRO')
		plt.title( 'Eps = ' + str(eps_est[n]) + '; MAE = ' + str(base_MAE[n]))
		plt.show()

	plt.plot(trace, Z_surf)
	plt.plot(trace, profs[best_ind,:].flatten())
	plt.plot(trace, Z_surf - 0.3 / np.sqrt(eps_est[best_ind]) * delta_TWT/2, 'k--')
	plt.plot(trace, Z_surf - delta_TWT *0.3/2)
	plt.plot(plains_trace, plains_z)
	plt.legend(['Surf','Generated Base','Depth-Corrected Reflector','Uncorrected Reflector','Plains'])
	plt.xlabel('Trace')
	plt.ylabel('Meters below MRO')
	plt.title('Best Eps = ' + str(eps_est[best_ind]) + '; MAE = ' + str(base_MAE[best_ind]))
	plt.show()

	plt.plot(eps_est, base_MAE, 'k.')
	plt.xlabel('Epsilon')
	plt.ylabel('MAE')
	plt.show()
	
	return

if __name__ == '__main__':
	main() 
