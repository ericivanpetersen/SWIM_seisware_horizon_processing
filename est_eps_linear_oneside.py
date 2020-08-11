import sys
import os
from mars_projections import *
from SWIM_horizons import *

# Path to MOLA geotiff file; must be changed if using script on new or different computer:
MOLA_file='/Users/eric/Documents/orig/supl/MOLA/DEM_global_mola128ppd_merged_mola64ppd/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif'

def main():
	# Usage: python ./est_eps_ATM.py {horizon file} {region} {sub_horiz} {surf_horiz} {Optional Arg: Orbit List}

	# Test if MOLA File exists, as hardwired into the script:
	file_exists = os.path.isfile(MOLA_file)
	if file_exists:
		print("MOLA File Found")
	else:
		print("ERROR: MOLA file not found; edit MOLA_file variable in script")
		return 

	#Test for correct number of inputs:
	if len(sys.argv) < 5:
		print()
		print('ERROR: Incorrect number of args (4 expected)')
		print()
		print('Usage: "python ./est_eps_linear_oneside.py {horizon file} {region} {sub_horiz} {surf_horiz}')
		print('horizon file = Seisware Export')
		print('region = onilus, utopia, arcadia, acidalia, or box[1-7]')
		print('sub_horiz = subsurface horizon')
		print('surf_horiz = surface horizon')
		print()
		sys.exit(1)

	horiz_file = sys.argv[1]
	region = sys.argv[2]
	sub_horiz = sys.argv[3]
	surf_horiz = sys.argv[4]
	if len(sys.argv) > 5:
		orbit_list = np.array([sys.argv[nn] for nn in range(5,len(sys.argv),1)], dtype='int')
	else:
		orbit_list = []
	
	if region == 'onilus': region = onilus
	if region == 'acidalia': region = acidalia
	if region == 'utopia': region = utopia
	if region == 'arcadia': region = arcadia
	if region == 'e_hellas': region = e_hellas
	if region == 'box2': region = box2
	if region == 'sirenum': region = sirenum
	if region == 'chryse': region = chryse
	if region == 'tharsis_tempe': region = tharsis_tempe
	if region == 'w_hellas': region = w_hellas
	if region == 'cimmeria': region = cimmeria

	data = swim_horizons(horiz_file, region)

	file_path = os.path.dirname(horiz_file) + '/Dielectric_Estimate/'

	data.estimate_epsilon_along_track_MOLA_minima(sub_horiz, surf_horiz, MOLA_file, file_path, orbit_list)

if __name__ == '__main__':
	main()
