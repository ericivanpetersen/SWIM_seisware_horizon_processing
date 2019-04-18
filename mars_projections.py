import numpy as np
from osgeo import osr

# Define Projections:

## Albers Equal-Area Projections Defined for SWIM Study regions in Seisware Projects:
onilus = 'PROJCS["Mars_Albers_AUTO", GEOGCS["Mars 2000", DATUM["D_Mars_2000",SPHEROID["Mars_2000_IAU_IAG",3396000.0,169.8]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",35],PARAMETER["Standard_Parallel_1",27.5],PARAMETER["Standard_Parallel_2",57.5],PARAMETER["Latitude_Of_Center",45],PARAMETER["Longitude_Of_Center",35],UNIT["Meter",1]]'
acidalia = 'PROJCS["Mars_Albers_AUTO", GEOGCS["Mars 2000", DATUM["D_Mars_2000",SPHEROID["Mars_2000_IAU_IAG",3396000.0,169.8]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",-35],PARAMETER["Standard_Parallel_1",27.5],PARAMETER["Standard_Parallel_2",57.5],PARAMETER["Latitude_Of_Center",45],PARAMETER["Longitude_Of_Center",-35],UNIT["Meter",1]]'
arcadia = 'PROJCS["Mars_Albers_AUTO", GEOGCS["Mars 2000", DATUM["D_Mars_2000",SPHEROID["Mars_2000_IAU_IAG",3396000.0,169.8]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",-172.5],PARAMETER["Standard_Parallel_1",27.5],PARAMETER["Standard_Parallel_2",57.5],PARAMETER["Latitude_Of_Center",45],PARAMETER["Longitude_Of_Center",-172.5],UNIT["Meter",1]]'
utopia = 'PROJCS["Mars_Albers_AUTO", GEOGCS["Mars 2000", DATUM["D_Mars_2000",SPHEROID["Mars_2000_IAU_IAG",3396000.0,169.8]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",110],PARAMETER["Standard_Parallel_1",27.5],PARAMETER["Standard_Parallel_2",57.5],PARAMETER["Latitude_Of_Center",45],PARAMETER["Longitude_Of_Center",110],UNIT["Meter",1]]'

# The Classic Mars 2000 Geographic Coordinate System:
mars_2000 = 'GEOGCS["Mars 2000",DATUM["D_Mars_2000",SPHEROID["Mars_2000_IAU_IAG",3396190.0,169.89444722361179]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]]'

# The Class Mars Equidistant Cylindrical Projection:
mars_equidistant_cylindrical = '+proj=eqc +lat_0=0 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=3396190 +b=3376200 +units=m +no_defs'

def define_proj(proj_code):
	"""Imports projection information for 
	projection code, whether it is EPSG, WKT,
	or Proj4 format.

	:param proj_code: projection code input"""

	# Create Projection:
	proj_srs = osr.SpatialReference()
	# Determine if EPSG Code, Well-Known Text (WKT), or Proj4:
	try:
		proj_num = int(proj_code)
		proj_srs.ImportFromEPSG(proj_num) #Try EPSG Code
	except ValueError:
		fs = proj_code[0]
		if(fs == 'G' or fs == 'P'):
			proj_srs.ImportFromWkt(proj_code) #WKT
		else:
			proj_srs.ImportFromProj4(proj_code) #Proj4

	return proj_srs

def create_transform(src, dst):
	"""Creates an OSR transformation.
	Input geometry definitions can be defined by:
		EPSG Codes
		Well-Known Text (WKT)
		Proj4

	:param src: Source geometry
	:param dst: Destination geometry
	:return: osr.CoordinateTransformation
	"""
	# Define Source Projection:
	src_srs = define_proj(src)

	# Define Target/Destination Projection:
	dst_srs = define_proj(dst)

	# Define the Transformation:
	transformation = osr.CoordinateTransformation(src_srs, dst_srs)

	return transformation

def transform_points(points, src, dst=mars_2000):
	"""Transform the coordinate reference system of a
	list of coordinates (points).
	Input geometries can be in EPSG Code, WKT, or Proj4
	
	:param float points: points to transform
	:param str src: Source Geometry
	:param str dst: Destination Geometry
	"""
	transform = create_transform(src, dst)
	points = transform.TransformPoints(points)
	return points

def transform_XY(x, y, src, dst=mars_2000):
	"""Same as transform points, but takes individual
	X & Y arrays as inputs.

	:param float x: x array
	:param float y: y array
	:param str src: source geometry
	:param str dst: destination geometry
	:output float x2: x coordinates in new geometry
	:output float y2: y coordinates in new geometry
	"""
	xy = []
	x2 = []
	y2 = []
	for nn in range(len(x)):
		xy.append([x[nn], y[nn]])
	xy2 = transform_points(xy, src, dst)
	for nn in range(len(x)):
		point = xy2[nn]
		x2.append(point[0])
		y2.append(point[1])
	x2 = np.array(x2)
	y2 = np.array(y2)

	return x2, y2
