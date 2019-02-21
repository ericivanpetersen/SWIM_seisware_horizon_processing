# SWIM_seisware_horizon_processing

Scripts and functions written in python that read in Seisware horizon exports, manage Mars-specific projections, process and produce data products relevent to the SWIM study.

## Dependencies

numpy, gdal, osgeo, os, sys, csv

## How to

To read in Seisware horizons and export csv files for each, use the following command:
```
"python make_refl_maps.py {horizon_export_file} {region}"
```
Where region = onilus, utopia, arcadia, or acidalia.



If you desire to make a reflection confidence map (0 or 1), use the following command:
```
"python make_refl_maps.py {horizon_export_file} {region} {subsurface horizon 1} {subsurface horizon 2} (etc.)"
```
for n subsurface horizons that you want to be represented as a "1" value in the map.

## More info:

mars_projections.py: useful projection information and functions for handling Mars projections, particularly those related to SWIM

seisware_horizons.py: defines how to read seisware horizons exports and define a class containing horizon information

make_refl_maps.py: customizes seisware_horizons.py for the SWIM project
