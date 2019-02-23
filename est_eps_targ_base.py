import sys
import os
from mars_projections import *
from SWIM_horizons import *

def main():
	# Usage: python ./est_eps_targ_base.py {horizon file} {region} {sub_horiz} {surf_horiz} {targ_horiz}

	#Test for correct number of inputs:
	if len(sys.argv) != 6:
		print
		print 'ERROR: Incorrect number of args (5 expected)'
		print 
		print 'Usage: "python ./est_eps_targ_base.py {horizon file} {region} {sub_horiz} {surf_horiz} {targ_horiz}'
		print 'horizon file = Seisware Export'
		print 'region = onilus, utopia, arcadia, or acidalia'
		print 'sub_horiz = subsurface horizon'
		print 'surf_horiz = surface horizon'
		print 'targ_horiz = target horizon'
		print
		sys.exit(1)

	horiz_file = sys.argv[1]
	region = sys.argv[2]
	sub_horiz = sys.argv[3]
	surf_horiz = sys.argv[4]
	targ_horiz = sys.argv[5]

	if region == 'onilus': region = onilus
	if region == 'acidalia': region = acidalia
        if region == 'utopia': region = utopia
        if region == 'arcadia': region = arcadia

	data = swim_horizons(horiz_file, region)

	pick_path = os.path.dirname(horiz_file) + '/'

	data.estimate_epsilon_targ_horiz(sub_horiz, surf_horiz, targ_horiz, pick_path)

if __name__ == '__main__':
	main()
