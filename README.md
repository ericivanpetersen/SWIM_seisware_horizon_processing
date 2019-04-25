# SWIM_seisware_horizon_processing

Scripts and functions written in python that read in Seisware horizon exports, manage Mars-specific projections, process and produce data products relevent to the SWIM study.

## Dependencies

numpy, gdal, osgeo, os, sys, csv, matplotlib

## How to

To estimate the real dielectric constant using the "Along-Track MOLA Minima Extrapolation Method (ATM)," obtain a geotiff of Mars MOLA elevation in an equidistant cylindrical projection, modify MOLA_file on line 7 in "est_eps_ATM.py" to direct the script to that geotiff filepath, and run the following command:
```
python est_eps_ATM.py {horizon_export_file} {region} {subsurface horizon} {surface horizon} {Optional additional arguments: orbit numbers to examine}
```
Without the additional arguments the code will loop you through all of the orbits containing the subsurface horizon. With the additional arguments it will only loop through the provided orbit numbers. In either case the code will save results in several files in a folder "/Dielectric_ATM/" in the folder wheere the {horizon_export_file} is provided. These include a "_results.csv" file containing trace-by-trace results, and a "_summary.csv" file providing results summarized for segments of subsurface reflector. When the code is completed it will also save a “_depth_result.csv” file which records the depths at each trace calculated using the median epsilon derived for the whole dataset.

When the code is running it will provide a MOLA profile of the track in question for you to scrutinize. You will then be asked to provide the number of segments between MOLA minima you want to find (0 indicates passing on the profile), then the range of latitudes to search for those minima. After defining the segment minima you will have a chance to review your work and decide to continue, keeping the segment or discarding, or retry with new defined minima.

If the code fails partway through your work, the “_results.csv” and “_summary.csv” files will be saved up until that point; you can restart the script and tell it which orbit you left off at.


Example 1 (Do all SHARAD Orbits):
```
python ./est_eps_ATM.py ../Horizon_Export/2019_04_08.txt onilus plains_sub1_EP surf_EP
```
Example 2 (Do only these two SHARAD Orbits, 187301 and 210401):
```
python ./est_eps_ATM.py ../Horizon_Export/2019_04_08.txt onilus plains_sub1_EP surf_EP 187301 210401
```
## More info:

mars_projections.py: useful projection information and functions for handling Mars projections, particularly those related to SWIM

general_functions.py: functions to aid in code abstraction for the purposes of writing files, radar and SWIM-related calculations

seisware_horizons.py: defines how to read seisware horizons exports and define a class containing horizon information

SWIM_horizons.py: customizes "seisware_horizons.py" for the SWIM project, particularly by the addition of dielectric estimation and derivation of consistency values.

