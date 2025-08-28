import rasterio as rio
import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import shapefile
import datetime
from pathlib import Path
from shapely.geometry import Point, shape
import xml.etree.ElementTree as ET
from urllib.request import urlretrieve

############################################
# Start of variables to change for your area
############################################

# Main GMRT (or other WGS84 digital elevation map in GeoTIFF format) file to load
# If using elevation data that is not from GMRT, the readme.txt section at the end
# should be updated to credit the correct source.
gmrt_file = 'GMRTv4_3_1_20250814topo.tif'

# Required output resolution. Ideally square, and (2^n + 1)
output_samples_long = 1025
output_samples_lat = 1025

# Location to output to
output_folder = Path('Output_PortsmouthHarbour')

################################
# General configuration settings
################################

# Do we want to replace nan and inf values (required for 5.10.? and earlier)
replace_nans = True

# If we want to get navigation aid (e.g. buoy) information from OSM
use_osm_map = True

# If we want to add coastline detail with the OSM coastline data
# Warning, this is currently very slow for reasonable size output samples
use_osm_coastline = True

# The OSM coastlines (land polygons) file to use
# This can be downloaded from
# https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip
osm_land_file = 'land-polygons-split-4326.zip'

# Values to use when applying overlap from OpenStreetMap coastlines file
min_sea_depth = 10
min_land_height = 1

# If base data is higher than filter_land_height_limit, or deeper than
# filter_sea_depth_limit, then don't check against OSM data (for speed)
filter_sea_depth_limit = 10
filter_land_height_limit = 10

#######################################
# End of settings, start of main script
#######################################

print('Script start: ' + str(datetime.datetime.now()))

# If output folder does not exist, make it
if not output_folder.exists():
    output_folder.mkdir()

# Load the data and show basic information
gmrt_tiff = rio.open(gmrt_file)
print("Data bounds: " + str(gmrt_tiff.bounds))
print("Data width: " + str(gmrt_tiff.width))
print("Data height: " + str(gmrt_tiff.height))
print("Output pixels per input (H): " + str(output_samples_long/gmrt_tiff.width))
print("Output pixels per input (V): " + str(output_samples_lat/gmrt_tiff.height))

# Get the required parameters
terrain_long = gmrt_tiff.bounds.left
terrain_lat = gmrt_tiff.bounds.bottom
terrain_long_extent = gmrt_tiff.bounds.right - gmrt_tiff.bounds.left
terrain_lat_extent = gmrt_tiff.bounds.top - gmrt_tiff.bounds.bottom

# Sample the data as required and load
out_shape = (output_samples_lat,output_samples_long)
resampling = rio.enums.Resampling.bilinear
gmrt_array = gmrt_tiff.read(indexes=1, out_shape=out_shape, resampling=resampling) # Note that this can read a different window and resample if required

# Replace nan and inf values
if replace_nans:
    gmrt_array = np.nan_to_num(gmrt_array, nan=-999, posinf=999, neginf=-999)

# Make a copy (instead of a view)
gmrt_array_working = gmrt_array.copy()

# Find the max sea depth and land max height
sea_max_depth = -1 * np.min(gmrt_array_working, where=~np.isnan(gmrt_array_working), initial=0)
terrain_max_height = np.max(gmrt_array_working, where=~np.isnan(gmrt_array_working), initial=0)

# Create terrain.ini file
terrain_ini_file_name = output_folder / 'terrain.ini'
terrain_ini_file = open(terrain_ini_file_name, 'w')

# General information - for now only one terrain map
terrain_ini_file.write("Number=1\n")
terrain_ini_file.write("MapImage=map.png\n")

terrain_ini_file.write("HeightMap(1)=height.f32\n")
terrain_ini_file.write("Texture(1)=texture.png\n")
terrain_ini_file.write("TerrainLong(1)=" + str(terrain_long) + "\n")
terrain_ini_file.write("TerrainLat(1)=" + str(terrain_lat) + "\n")
terrain_ini_file.write("TerrainLongExtent(1)=" + str(terrain_long_extent) + "\n")
terrain_ini_file.write("TerrainLatExtent(1)=" + str(terrain_lat_extent) + "\n")
# These are required for the .f32 format
terrain_ini_file.write("TerrainHeightMapRows(1)=" + str(output_samples_lat) + "\n")
terrain_ini_file.write("TerrainHeightMapColumns(1)=" + str(output_samples_long) + "\n")
# The details below will actually be read from the .f32 file, they are just here for information:
terrain_ini_file.write("TerrainMaxHeight(1)=" + str(terrain_max_height) + "\n")
terrain_ini_file.write("SeaMaxDepth(1)=" + str(sea_max_depth) + "\n")
terrain_ini_file.close()

if use_osm_coastline:
    # Find the bounding box of the area of interest
    world_bounding_box = (terrain_long,
                          terrain_lat,
                          terrain_long + terrain_long_extent,
                          terrain_lat + terrain_lat_extent)

    # Load the OpenStreetMap land polygon shapefile
    print('Loading OpenStreetMap land polygon shapefile\n')
    shp = shapefile.Reader(osm_land_file) #open the shapefile
    filtered_shapes = shp.shapes(bbox=world_bounding_box) # get all the polygons in our area
    number_of_shapes = len(filtered_shapes)

    filter_sample = 0
    total_samples = output_samples_long * output_samples_lat
    for x_index in range(output_samples_long):
        for y_index in range(output_samples_lat):
            if filter_sample % 1000 == 0:
                completion_percent = 100 * filter_sample / total_samples
                print(f'Filtering points {completion_percent:.2f}%')
            base_data_height = gmrt_array_working[y_index, x_index]
            if (base_data_height < filter_land_height_limit) and (base_data_height > -1 * filter_sea_depth_limit):
                # Heights in range to check against OSM data
                x_sample = terrain_long + x_index * (terrain_long_extent / output_samples_long)
                y_sample = terrain_lat + terrain_lat_extent - y_index * (terrain_lat_extent / output_samples_lat)
                # TODO: Check if these are off by 1 at max
                point_to_check = (x_sample, y_sample)
                point_is_land = False

                for shape_index in range(number_of_shapes):
                    if Point(point_to_check).within(shape(filtered_shapes[shape_index])):
                        point_is_land = True
                        # No need to check more polygons
                        break
                if point_is_land:
                    gmrt_array_working[y_index, x_index] = max(base_data_height, min_land_height)
                else:
                    gmrt_array_working[y_index, x_index] = min(base_data_height, -1 * min_sea_depth)
            filter_sample = filter_sample + 1

# Create the .f32 file
terrain_f32_file_name = output_folder / 'height.f32'
terrain_f32_file = open(terrain_f32_file_name, 'wb')
gmrt_array_flipped = np.flipud(gmrt_array_working)
gmrt_array_flipped.tofile(terrain_f32_file)
terrain_f32_file.close()

# Create a map file
terrain_map_file_name = output_folder / 'map.png'
img.imsave(terrain_map_file_name, gmrt_array_working, format='png', vmin=-10, vmax=10)

# Create a texture file
terrain_texture_file_name = output_folder / 'texture.png'
img.imsave(terrain_texture_file_name, gmrt_array_working, cmap=plt.get_cmap('summer'), format='png', vmin=0, vmax=terrain_max_height)

# End of main map generation

# Start of processing OSM map file for buoys, lights etc

# Utility function to get value from child if present
def get_child_value(xml_item, key):
    result = ''
    for xml_child in xml_item:
        if 'k' in xml_child.keys() and 'v' in xml_child.keys() and xml_child.attrib['k'] == key:
            result = xml_child.attrib['v']
    return result

if use_osm_map:
    print('Downloading OpenStreetMap data for area\n')
    osm_map_file = 'map.osm'
    osm_map_url = ('https://overpass-api.de/api/map?bbox=' +
                   str(terrain_long) + ',' +
                   str(terrain_lat) + ',' +
                   str(terrain_long + terrain_long_extent) + ',' +
                   str(terrain_lat + terrain_lat_extent))

    # Download the OSM Map here (TODO: Make this optional)
    urlretrieve(osm_map_url, osm_map_file)

    # create element tree object
    tree = ET.parse(osm_map_file)
    # get root element
    root = tree.getroot()

    buoys = []

    for item in root.findall('node'):
        if len(item) > 0:
            seamark_type = get_child_value(item, 'seamark:type')
            if seamark_type != '':
                seamark_category = get_child_value(item, 'seamark:' + seamark_type + ':category')
                seamark_light_character = get_child_value(item, 'seamark:light:character')
                seamark_light_colour = get_child_value(item, 'seamark:light:colour')
                seamark_lat = item.attrib['lat']
                seamark_lon = item.attrib['lon']
                # Also to look at - sequence, period, height, range
                # Also look at sector lights, which will need to be handled differently

                # Add known buoy types for now
                if seamark_type == 'beacon_lateral' and seamark_category == 'port':
                    buoys.append(['port_post', seamark_lon, seamark_lat, True])
                if seamark_type == 'beacon_lateral' and seamark_category == 'starboard':
                    buoys.append(['stbd_post', seamark_lon, seamark_lat, True])

                if seamark_type == 'beacon_special_purpose':
                    buoys.append(['special_post', seamark_lon, seamark_lat, True])

                if seamark_type == 'buoy_special_purpose' and seamark_category == 'mooring':
                    buoys.append(['mooring', seamark_lon, seamark_lat, False])

                if seamark_type == 'buoy_lateral' and seamark_category == 'port':
                    buoys.append(['port', seamark_lon, seamark_lat, False])
                if seamark_type == 'buoy_lateral' and seamark_category == 'starboard':
                    buoys.append(['stbd', seamark_lon, seamark_lat, False])

                if seamark_type == 'buoy_cardinal' and seamark_category == 'north':
                    buoys.append(['north', seamark_lon, seamark_lat, False])
                if seamark_type == 'buoy_cardinal' and seamark_category == 'east':
                    buoys.append(['east', seamark_lon, seamark_lat, False])
                if seamark_type == 'buoy_cardinal' and seamark_category == 'south':
                    buoys.append(['south', seamark_lon, seamark_lat, False])
                if seamark_type == 'buoy_cardinal' and seamark_category == 'west':
                    buoys.append(['west', seamark_lon, seamark_lat, False])

                # TODO: beacon_cardinal needed as well, plus all other buoy types

    # Create buoy.ini file
    buoy_ini_file_name = output_folder / 'buoy.ini'
    buoy_ini_file = open(buoy_ini_file_name, 'w')

    buoy_ini_file.write('Number=' + str(len(buoys)) + '\n\n')
    for i, buoy in enumerate(buoys):
        buoy_ini_file.write('Type(' + str(i+1) + ')=' + buoy[0] + '\n' )
        buoy_ini_file.write('Long(' + str(i + 1) + ')=' + buoy[1] + '\n')
        buoy_ini_file.write('Lat(' + str(i + 1) + ')=' + buoy[2] + '\n')
        if buoy[3]:
            buoy_ini_file.write('Grounded(' + str(i + 1) + ')=1\n')
        buoy_ini_file.write('\n')

    buoy_ini_file.close()

# Write a readme.txt file (For GMRT and OSM data)
readme_file_name = output_folder / 'readme.txt'
readme_file = open(readme_file_name, 'w')
readme_file.write('Elevation and Bathymetry data from the Global Multi-Resolution Topography Synthesis (GMRT)\n')
readme_file.write('For details, please see Ryan, W. B. F., S.M. Carbotte, J. Coplan, S. O\'Hara, A. Melkonian, R. Arko,\n'
                  'R.A. Weissel, V. Ferrini, A. Goodwillie, F. Nitsche, J. Bonczkowski, and R. Zemsky (2009),\n'
                  'Global Multi-Resolution Topography (GMRT) synthesis data set, Geochem. Geophys. Geosyst., \n'
                  '10, Q03014, doi:10.1029/2008GC002332.\n'
                  '\n'
                  'Coastline data and navigation aid data from from OpenStreetMap data, available\n'
                  'under the Open Database License.\n'
                  '\n'
                  'NOT FOR USE IN NAVIGATION\n')
readme_file.close()

# End
print('Script end: ' + str(datetime.datetime.now()))
