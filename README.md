# SWIM_seisware_horizon_processing

Scripts and functions written in python that read in Seisware horizon exports, manage Mars-specific projections, process and produce data products relevent to the SWIM study.

## Dependencies

numpy, gdal, osgeo, os, sys, csv

## How to

To read in Seisware horizons and export csv files for each, use the following command:
```
python make_refl_maps.py {horizon_export_file} {region}
```
Where region = onilus, utopia, arcadia, or acidalia.



If you desire to make a reflection confidence map (0 or 1), use the following command:
```
python make_refl_maps.py {horizon_export_file} {region} {subsurface horizon 1} {subsurface horizon 2} (etc.)
```
for n subsurface horizons that you want to be represented as a "1" value in the map.


To estimate the dielectric constant by depth-correcting a subsurface reflector to a target horizon, use the following command:
```
python est_eps_targ_base.py {horizon_export_file} {region} {subsurface horizon} {surface horizon} {target horizon}
```

To estimate the real dielectric constant using the "Along-Track MOLA Minima Extrapolation Method (ATM)," obtain a geotiff of Mars MOLA elevation in an equidistant cylindrical projection, modify MOLA_file on line 7 in "est_eps_ATM.py" to direct the script to that geotiff filepath, and run the following command:
```
python est_eps_ATM.py {horizon_export_file} {region} {subsurface horizon} {surface horizon} {Optional additional arguments: orbit numbers to examine}
```
Without the additional arguments the code will loop you through all of the orbits containing the subsurface horizon. With the additional arguments it will only loop through the provided orbit numbers. In either case the code will save results in several files in a folder "/Dielectric_ATM/" in the folder wheere the {horizon_export_file} is provided. These include a "_result.csv" file containing trace-by-trace results, and a "_summary.csv" file providing results summarized for segments of subsurface reflector.
When the code is running it will provide a MOLA profile of the track in question for you to scrutinize. You will then be asked to provide the number of segments between MOLA minima you want to find (0 indicates passing on the profile), then the range of latitudes to search for those minima. After defining the segment minima you will have a chance to review your work and decide to continue, keeping the segment or discarding, or retry with new defined minima.
In its current form this code is not conducive to doing the work in bite-sized chunks; need to run the whole thing to the end to save the trace-by-trace "_results.csv," for example. Thus I advise talking notes on your inputs as you work through it. And if you do it in chunks using the optional orbit # arguments, be sure to change the output filenames so you don't overwrite them.

## More info:

mars_projections.py: useful projection information and functions for handling Mars projections, particularly those related to SWIM

seisware_horizons.py: defines how to read seisware horizons exports and define a class containing horizon information

SWIM_horizons.py: customizes "seisware_horizons.py" for the SWIM project, particularly by the addition of the consistency values

make_refl_maps.py: can make reflector consistency maps; also produces simple .csv outputs of individual horizons
