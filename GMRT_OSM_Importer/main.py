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

# Default light settings
default_light_height = 5  # m
default_light_range = 5  # Nm

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
    osm_map_file = 'map.osm'
    osm_map_url = ('https://overpass-api.de/api/map?bbox=' +
                   str(terrain_long) + ',' +
                   str(terrain_lat) + ',' +
                   str(terrain_long + terrain_long_extent) + ',' +
                   str(terrain_lat + terrain_lat_extent))

    # Download the OSM Map here (TODO: Make this optional, and allow force override)
    if not Path.is_file(Path(osm_map_file)):
        print('Downloading OpenStreetMap data for area from ' + osm_map_url + '\n')
        urlretrieve(osm_map_url, osm_map_file)
    else:
        print('Using cached ' + osm_map_file + ' for OpenStreetMap data for area\n' )

    # create element tree object
    tree = ET.parse(osm_map_file)
    # get root element
    root = tree.getroot()

    # Iterate for buoys/posts and similar
    buoys = []
    for node_index, item in enumerate(root.findall('node')):
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
                    buoys.append(['port_post', seamark_lon, seamark_lat, True, node_index])
                if seamark_type == 'beacon_lateral' and seamark_category == 'starboard':
                    buoys.append(['stbd_post', seamark_lon, seamark_lat, True, node_index])

                if seamark_type == 'beacon_special_purpose':
                    buoys.append(['special_post', seamark_lon, seamark_lat, True, node_index])

                if seamark_type == 'buoy_special_purpose' and seamark_category == 'mooring':
                    buoys.append(['mooring', seamark_lon, seamark_lat, False, node_index])

                if seamark_type == 'buoy_lateral' and seamark_category == 'port':
                    buoys.append(['port', seamark_lon, seamark_lat, False, node_index])
                if seamark_type == 'buoy_lateral' and seamark_category == 'starboard':
                    buoys.append(['stbd', seamark_lon, seamark_lat, False, node_index])

                if seamark_type == 'buoy_cardinal' and seamark_category == 'north':
                    buoys.append(['north', seamark_lon, seamark_lat, False, node_index])
                if seamark_type == 'buoy_cardinal' and seamark_category == 'east':
                    buoys.append(['east', seamark_lon, seamark_lat, False, node_index])
                if seamark_type == 'buoy_cardinal' and seamark_category == 'south':
                    buoys.append(['south', seamark_lon, seamark_lat, False, node_index])
                if seamark_type == 'buoy_cardinal' and seamark_category == 'west':
                    buoys.append(['west', seamark_lon, seamark_lat, False, node_index])

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

    # Iterate for lights
    lights = []
    for node_index, item in enumerate(root.findall('node')):
        if len(item) > 0:
            # Bodge to find potential multiple (e.g. sector) lights
            for multi in range(10):
                if multi == 0:
                    multi_str = ''
                else:
                    multi_str = str(multi) + ':'

                light_character = get_child_value(item, 'seamark:light:' + multi_str + 'character')
                if light_character != '':
                    # Some type of light found

                    # Find if there's a matching entry in the 'buoys' list (same node_index)
                    light_buoy = 0  # Initially assume not found
                    for buoy_index, buoy in enumerate(buoys):
                        if buoy[4] == node_index:
                            light_buoy = buoy_index + 1

                    light_dict = {
                        "light_lat" : item.attrib['lat'],
                        "light_lon" : item.attrib['lon'],
                        "light_colour" : get_child_value(item, 'seamark:light:' + multi_str + 'colour'),
                        "light_range" : get_child_value(item, 'seamark:light:' + multi_str + 'range'),
                        "light_period" : get_child_value(item, 'seamark:light:' + multi_str + 'period'),
                        "light_group" : get_child_value(item, 'seamark:light:' + multi_str + 'group'),
                        "light_character" : light_character,
                        "sector_start" : get_child_value(item, 'seamark:light:' + multi_str + 'sector_start'),
                        "sector_end": get_child_value(item, 'seamark:light:' + multi_str + 'sector_end'),
                        "light_height": get_child_value(item, 'seamark:light:' + multi_str + 'height'),
                        "light_buoy": light_buoy,
                        "node_index" : node_index
                    }
                    lights.append(light_dict)

    # Create light.ini file
    light_ini_file_name = output_folder / 'light.ini'
    light_ini_file = open(light_ini_file_name, 'w')

    light_sequence = {
        'F' : 'LLLL',
        'Fl' : 'LLDDDDDD',
        'Q' : 'LLDDD',
        'VQ' : 'LDD',
        'LFl' : 'LLLLLLLLDDDDDDDDDDDD',
        'Iso' : 'LLLDDD',
        'Oc' : 'LLLLLLDD'
    }

    light_ini_file.write('Number=' + str(len(lights)) + '\n\n')
    for i, light in enumerate(lights):

        if light["light_character"] in light_sequence:
            sequence = light_sequence[light["light_character"]]

            # Duplicate for group
            if light["light_group"] != '':
                if int(light["light_group"]) > 1:
                    sequence = sequence * int(light["light_group"])

            # Pad to length
            if light["light_period"] != '':
                while len(sequence) * 0.25 < float(light["light_period"]):
                    if light["light_character"] == 'Oc':
                        # Special case for occulting, pad with light
                        sequence = sequence + "L"
                    else:
                        # Normal case
                        sequence = sequence + "D"

            if light["light_range"] != '':
                light_range = light["light_range"]
            else:
                # Default value, nautical miles
                light_range = str(default_light_range)

            light_ini_file.write('Desc(' + str(i + 1) + ')=' + str(light) + '\n')

            # Define position, either through lat/long, or which buoy/mark it's associated with
            if light["light_buoy"] > 0:
                light_ini_file.write('Buoy(' + str(i + 1) + ')=' + str(light["light_buoy"]) + '\n')
            else:
                light_ini_file.write('Long(' + str(i + 1) + ')=' + light["light_lon"] + '\n')
                light_ini_file.write('Lat(' + str(i + 1) + ')=' + light["light_lat"] + '\n')

            light_ini_file.write('Sequence(' + str(i + 1) + ')=' + sequence + '\n')
            light_ini_file.write('Range(' + str(i + 1) + ')=' + light_range + '\n')

            # For phase, use the node_index a bit like a hash input,
            # so lights on same node with the same length sequence are in phase
            phase_start = (light["node_index"] % len(sequence)) + 1
            light_ini_file.write('PhaseStart(' + str(i + 1) + ')=' + str(phase_start) + '\n')

            if light["light_height"] != '':
                # Note this is actually height above mean sea level from OSM definition,
                # might need adjusting to be above chart datum
                light_ini_file.write(
                    'Height(' + str(i + 1) + ')=' + light["light_height"] + '\n')
                light_ini_file.write('Absolute(' + str(i + 1) + ')=1\n')
            else:
                light_ini_file.write(
                    'Height(' + str(i + 1) + ')=' + str(default_light_height) + '\n')

            if light["light_colour"] == 'white':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '255' + '\n')
            elif light["light_colour"] == 'red':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '0' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '0' + '\n')
            elif light["light_colour"] == 'green':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '0' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '0' + '\n')
            elif light["light_colour"] == 'blue':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '0' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '0' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '255' + '\n')
            elif light["light_colour"] == 'yellow':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '0' + '\n')
            elif light["light_colour"] == 'amber':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '128' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '0' + '\n')
            elif light["light_colour"] == 'orange': # Currently same as amber
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '128' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '0' + '\n')
            elif light["light_colour"] == 'violet':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '128' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '255' + '\n')
            elif light["light_colour"] == 'magenta':
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '128' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '0' + '\n')
            else:
                print('Light colour ' + light["light_colour"] + ' not mapped\n')
                light_ini_file.write('Red(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Green(' + str(i + 1) + ')=' + '255' + '\n')
                light_ini_file.write('Blue(' + str(i + 1) + ')=' + '255' + '\n')

            if (light["sector_start"] != '') and (light["sector_end"] != ''):
                # 180 degree difference in definitions of light angles between OSM and Bridge Command (OSM is start bearing seen by observer)
                start_angle = float(light["sector_start"]) + 180
                end_angle = float(light["sector_end"]) + 180

                # Normalise angles
                while start_angle > 360:
                    start_angle = start_angle - 360
                    end_angle = end_angle - 360

                # Ensure end angle is always greater than start
                while end_angle < start_angle:
                    end_angle = end_angle + 360

                light_ini_file.write('StartAngle(' + str(i + 1) + ')=' + str(start_angle) + '\n')
                light_ini_file.write('EndAngle(' + str(i + 1) + ')=' + str(end_angle) + '\n')
            else:
                light_ini_file.write('StartAngle(' + str(i + 1) + ')=' + '0' + '\n')
                light_ini_file.write('EndAngle(' + str(i + 1) + ')=' + '360' + '\n')

            light_ini_file.write('\n')
        else:
            print('Light character ' + light["light_character"] + ' not mapped\n')

    light_ini_file.close()

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
