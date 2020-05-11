import os
import numpy as np
from SWIM_horizons import *
from mars_projections import *

### File paths; make sure these are correct for your computer:
MOLA_file='/Users/eric/Documents/orig/supl/MOLA/DEM_global_mola128ppd_merged_mola64ppd/mola128_mola64_merge_90Nto90S_SimpleC_clon0.tif' # MOLA File
datafile = '../../Horizon_Export/Box2/box2_2020_02_24.txt' # Data file with horizon export from Seisware
filepath = '../../Horizon_Export/Box2/' # Path to project folder (data files saved here)

def main():
    # Test if MOLA File exists, as hardwired into the script:
    file_exists = os.path.isfile(MOLA_file)
    if file_exists:
        print("MOLA File Found")
    else:
        print("ERROR: MOLA file not found; edit MOLA_file variable in script")
        return 

    orbit_selection = np.loadtxt('./calibration_orbit_selection.csv') # Read orbit selection list

    # Read in and manage horizons:
    data = swim_horizons(datafile, box2)
    data.combine_horizons('lda_surf_EP', 'plains_surf_EP', 'all_surf')

    # Perform dielectric estimation:
    data.estimate_epsilon_along_track_linear('lda_sub1_EP','all_surf',MOLA_file,filepath+'eps_lin_ext/',orbit_selection,method='polyfit')

    return

if __name__ == "__main__":
    main()