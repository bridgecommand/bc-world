import datetime
import os.path
import warnings
from pathlib import Path
import xml.etree.ElementTree as ET
from urllib.request import urlretrieve
import rasterio as rio
import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import shapefile
from shapely.geometry import Polygon, Point
from shapely.strtree import STRtree

############################################
# Start of variables to change for your area
############################################

# Main GMRT (or other WGS84 digital elevation map in GeoTIFF format) file to load
# If using elevation data that is not from GMRT, the readme.txt section at the end
# should be updated to credit the correct source.
gmrt_file = "GMRTv4_3_1_PortsmouthHarbour.tif"

IALA_region = "A"  # 'A' or 'B'

# Required output resolution. Ideally square, and (2^n + 1)
output_samples_long = 2049
output_samples_lat = 2049

# Location to output to, from GMRT filename, users can change this to a more readable path
output_folder = Path(gmrt_file + "_output")

################################
# General configuration settings
################################

# Do we want to replace nan and inf values (required for 5.10.? and earlier)
replace_nans = True

# If we want to get navigation aid (e.g. buoy) information from OSM
use_osm_map = True

# If we want to generate 3d models of buildings and similar on land from OSM (OSM2World)
# Note that this requires osmium and osm2world and java to be available
use_osm2world = True
# Divide world into grid to avoid excessively large models
x_divisions_OSM2World = 5
y_divisions_OSM2World = 5
osm_2_world_exe = Path("osm2world/OSM2World-latest-bin/osm2world-windows.bat")
java_home = Path("openjdk-25.0.1_windows-x64_bin/jdk-25.0.1")

# If we want to add coastline detail with the OSM coastline data
# Warning, this is currently very slow for reasonable size output samples
use_osm_coastline = True

# The OSM coastlines (land polygons) file to use
# This can be downloaded from
# https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip
osm_land_file = "land-polygons-split-4326.zip"
if use_osm_coastline and (not Path.is_file(Path(osm_land_file))):
    warnings.warn(
        "Open Street Map land polygons file not found! Continuing without coastline data.\nPlease download "
        "it from https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip"
    )
    use_osm_coastline = False

# Values to use when applying overlay from OpenStreetMap coastlines file
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

print("Script start: " + str(datetime.datetime.now()))

# If output folder does not exist, make it
if not output_folder.exists():
    output_folder.mkdir()

# Load the data and show basic information
gmrt_tiff = rio.open(gmrt_file)
print("Data bounds: " + str(gmrt_tiff.bounds))
print("Data width: " + str(gmrt_tiff.width))
print("Data height: " + str(gmrt_tiff.height))
print("Output pixels per input (H): " + str(output_samples_long / gmrt_tiff.width))
print("Output pixels per input (V): " + str(output_samples_lat / gmrt_tiff.height))

# Get the required parameters
terrain_long = gmrt_tiff.bounds.left
terrain_lat = gmrt_tiff.bounds.bottom
terrain_long_extent = gmrt_tiff.bounds.right - gmrt_tiff.bounds.left
terrain_lat_extent = gmrt_tiff.bounds.top - gmrt_tiff.bounds.bottom

# Sample the data as required and load
out_shape = (output_samples_lat, output_samples_long)
resampling = rio.enums.Resampling.bilinear
gmrt_array = gmrt_tiff.read(
    indexes=1, out_shape=out_shape, resampling=resampling
)  # Note that this can read a different window and resample if required

# Replace nan and inf values
if replace_nans:
    gmrt_array = np.nan_to_num(gmrt_array, nan=-999, posinf=999, neginf=-999)

# Make a copy (instead of a view)
gmrt_array_working = gmrt_array.copy()

# Find the max sea depth and land max height
sea_max_depth = -1 * np.min(
    gmrt_array_working, where=~np.isnan(gmrt_array_working), initial=0
)
terrain_max_height = np.max(
    gmrt_array_working, where=~np.isnan(gmrt_array_working), initial=0
)

# Create terrain.ini file
terrain_ini_file_name = output_folder / "terrain.ini"
terrain_ini_file = open(terrain_ini_file_name, "w")

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
    world_bounding_box = (
        terrain_long,
        terrain_lat,
        terrain_long + terrain_long_extent,
        terrain_lat + terrain_lat_extent,
    )

    # Load the OpenStreetMap land polygon shapefile
    print("Loading OpenStreetMap land polygon shapefile\n")
    shp = shapefile.Reader(osm_land_file)  # open the shapefile
    filtered_shapes = shp.shapes(
        bbox=world_bounding_box
    )  # get all the polygons in our area
    number_of_shapes = len(filtered_shapes)

    filter_sample = 0
    total_samples = output_samples_long * output_samples_lat
    samples = []
    for x_index in range(output_samples_long):
        for y_index in range(output_samples_lat):
            x_sample = terrain_long + x_index * (
                terrain_long_extent / (output_samples_long - 1)
            )
            y_sample = (
                terrain_lat
                + terrain_lat_extent
                - y_index * (terrain_lat_extent / (output_samples_lat - 1))
            )
            samples.append(Point(x_sample, y_sample))

    samples_tree = STRtree(samples)

    land_points = []
    for shape_index in range(number_of_shapes):
        land_poly = Polygon(filtered_shapes[shape_index].points)
        res = samples_tree.query(land_poly, predicate="contains")
        land_points.extend(samples_tree.geometries.take(res).tolist())

    is_land_array = np.zeros((output_samples_long, output_samples_lat))
    for land_point in land_points:
        land_point_x = land_point.x
        land_point_y = land_point.y
        land_point_index_x = round(
            (land_point_x - terrain_long)
            * (output_samples_long - 1)
            / terrain_long_extent
        )
        land_point_index_y = round(
            (terrain_lat + terrain_lat_extent - land_point_y)
            * (output_samples_lat - 1)
            / terrain_lat_extent
        )
        is_land_array[land_point_index_x, land_point_index_y] = 1

    for x_index in range(output_samples_long):
        for y_index in range(output_samples_lat):
            point_is_land = False
            base_data_height = gmrt_array_working[y_index, x_index]
            if is_land_array[x_index, y_index] > 0:
                point_is_land = True
            if point_is_land:
                gmrt_array_working[y_index, x_index] = max(
                    base_data_height, min_land_height
                )
            else:
                gmrt_array_working[y_index, x_index] = min(
                    base_data_height, -1 * min_sea_depth
                )

# Create the .f32 file
terrain_f32_file_name = output_folder / "height.f32"
terrain_f32_file = open(terrain_f32_file_name, "wb")
gmrt_array_flipped = np.flipud(gmrt_array_working)
gmrt_array_flipped.tofile(terrain_f32_file)
terrain_f32_file.close()

# Create a map file
terrain_map_file_name = output_folder / "map.png"
img.imsave(terrain_map_file_name, gmrt_array_working, format="png", vmin=-10, vmax=10)

# Create a texture file
terrain_texture_file_name = output_folder / "texture.png"
img.imsave(
    terrain_texture_file_name,
    gmrt_array_working,
    cmap=plt.get_cmap("summer"),
    format="png",
    vmin=0,
    vmax=terrain_max_height,
)

# End of main map generation

# Start of processing OSM map file for buoys, lights etc


# Utility function to get value from child if present
def get_child_value(xml_item, key):
    result = ""
    for xml_child in xml_item:
        if (
            "k" in xml_child.keys()
            and "v" in xml_child.keys()
            and xml_child.attrib["k"] == key
        ):
            result = xml_child.attrib["v"]
    return result


if use_osm_map:
    osm_map_file = gmrt_file + ".osm"
    osm_map_url = (
        "https://overpass-api.de/api/map?bbox="
        + str(terrain_long)
        + ","
        + str(terrain_lat)
        + ","
        + str(terrain_long + terrain_long_extent)
        + ","
        + str(terrain_lat + terrain_lat_extent)
    )

    # Download the OSM Map here
    if not Path.is_file(Path(osm_map_file)):
        print("Downloading OpenStreetMap data for area from " + osm_map_url)
        urlretrieve(osm_map_url, osm_map_file)
        print("Saved OpenStreetMap data to " + osm_map_file + "\n")
    else:
        print("Using cached " + osm_map_file + " for OpenStreetMap data for area.")

    # create element tree object
    tree = ET.parse(osm_map_file)
    # get root element
    root = tree.getroot()

    # Iterate for buoys/posts and similar
    buoys = []
    for node_index, item in enumerate(root.findall("node")):
        if len(item) > 0:
            seamark_type = get_child_value(item, "seamark:type")
            if seamark_type != "":
                seamark_category = get_child_value(
                    item, "seamark:" + seamark_type + ":category"
                )
                seamark_light_character = get_child_value(
                    item, "seamark:light:character"
                )
                seamark_light_colour = get_child_value(item, "seamark:light:colour")
                seamark_lat = item.attrib["lat"]
                seamark_lon = item.attrib["lon"]
                # Also to look at - sequence, period, height, range
                # Also look at sector lights, which will need to be handled differently

                # Add known buoy types for now
                if seamark_type == "beacon_lateral" and seamark_category == "port":
                    buoys.append(
                        ["port_post", seamark_lon, seamark_lat, True, node_index]
                    )
                if seamark_type == "beacon_lateral" and seamark_category == "starboard":
                    buoys.append(
                        ["stbd_post", seamark_lon, seamark_lat, True, node_index]
                    )

                if seamark_type == "beacon_special_purpose":
                    buoys.append(
                        ["special_post", seamark_lon, seamark_lat, True, node_index]
                    )

                if (
                    seamark_type == "buoy_special_purpose"
                    and seamark_category == "mooring"
                ):
                    buoys.append(
                        ["mooring", seamark_lon, seamark_lat, False, node_index]
                    )

                if seamark_type == "buoy_lateral" and seamark_category == "port":
                    buoys.append(
                        ["port_med", seamark_lon, seamark_lat, False, node_index]
                    )
                if seamark_type == "buoy_lateral" and seamark_category == "starboard":
                    buoys.append(
                        ["stbd_med", seamark_lon, seamark_lat, False, node_index]
                    )

                if seamark_type == "buoy_cardinal" and seamark_category == "north":
                    buoys.append(["north", seamark_lon, seamark_lat, False, node_index])
                if seamark_type == "buoy_cardinal" and seamark_category == "east":
                    buoys.append(["east", seamark_lon, seamark_lat, False, node_index])
                if seamark_type == "buoy_cardinal" and seamark_category == "south":
                    buoys.append(["south", seamark_lon, seamark_lat, False, node_index])
                if seamark_type == "buoy_cardinal" and seamark_category == "west":
                    buoys.append(["west", seamark_lon, seamark_lat, False, node_index])

                # TODO: beacon_cardinal needed as well, plus all other buoy types

    # Create buoy.ini file
    buoy_ini_file_name = output_folder / "buoy.ini"
    buoy_ini_file = open(buoy_ini_file_name, "w")

    # If IALA B, swap 'port' and 'stbd' for now:
    if IALA_region == "B":
        for buoy in buoys:
            if "stbd" in buoy[0]:
                buoy[0] = buoy[0].replace("stbd", "port")
            elif "port" in buoy[0]:
                buoy[0] = buoy[0].replace("port", "stbd")

    buoy_ini_file.write("Number=" + str(len(buoys)) + "\n\n")
    for i, buoy in enumerate(buoys):
        buoy_ini_file.write("Type(" + str(i + 1) + ")=" + buoy[0] + "\n")
        buoy_ini_file.write("Long(" + str(i + 1) + ")=" + buoy[1] + "\n")
        buoy_ini_file.write("Lat(" + str(i + 1) + ")=" + buoy[2] + "\n")
        if buoy[3]:
            buoy_ini_file.write("Grounded(" + str(i + 1) + ")=1\n")
        buoy_ini_file.write("\n")

    buoy_ini_file.close()

    # Iterate for lights
    lights = []
    for node_index, item in enumerate(root.findall("node")):
        if len(item) > 0:
            # Bodge to find potential multiple (e.g. sector) lights
            for multi in range(10):
                if multi == 0:
                    multi_str = ""
                else:
                    multi_str = str(multi) + ":"

                light_character = get_child_value(
                    item, "seamark:light:" + multi_str + "character"
                )
                if light_character != "":
                    # Some type of light found

                    # Find if there's a matching entry in the 'buoys' list (same node_index)
                    light_buoy = 0  # Initially assume not found
                    for buoy_index, buoy in enumerate(buoys):
                        if buoy[4] == node_index:
                            light_buoy = buoy_index + 1

                    light_dict = {
                        "light_lat": item.attrib["lat"],
                        "light_lon": item.attrib["lon"],
                        "light_colour": get_child_value(
                            item, "seamark:light:" + multi_str + "colour"
                        ),
                        "light_range": get_child_value(
                            item, "seamark:light:" + multi_str + "range"
                        ),
                        "light_period": get_child_value(
                            item, "seamark:light:" + multi_str + "period"
                        ),
                        "light_group": get_child_value(
                            item, "seamark:light:" + multi_str + "group"
                        ),
                        "light_character": light_character,
                        "sector_start": get_child_value(
                            item, "seamark:light:" + multi_str + "sector_start"
                        ),
                        "sector_end": get_child_value(
                            item, "seamark:light:" + multi_str + "sector_end"
                        ),
                        "light_height": get_child_value(
                            item, "seamark:light:" + multi_str + "height"
                        ),
                        "light_buoy": light_buoy,
                        "node_index": node_index,
                    }
                    lights.append(light_dict)

    # Create light.ini file
    light_ini_file_name = output_folder / "light.ini"
    light_ini_file = open(light_ini_file_name, "w")

    light_sequence = {
        "F": "LLLL",
        "Fl": "LLDDDDDD",
        "Q": "LLDDD",
        "VQ": "LDD",
        "LFl": "LLLLLLLLDDDDDDDDDDDD",
        "Iso": "LLLDDD",
        "Oc": "LLLLLLDD",
    }

    light_ini_file.write("Number=" + str(len(lights)) + "\n\n")
    for i, light in enumerate(lights):

        if light["light_character"] in light_sequence:
            sequence = light_sequence[light["light_character"]]

            # Duplicate for group
            if light["light_group"] != "":
                if int(light["light_group"]) > 1:
                    sequence = sequence * int(light["light_group"])

            # Pad to length
            if light["light_period"] != "":
                while len(sequence) * 0.25 < float(light["light_period"]):
                    if light["light_character"] == "Oc":
                        # Special case for occulting, pad with light
                        sequence = sequence + "L"
                    else:
                        # Normal case
                        sequence = sequence + "D"

            if light["light_range"] != "":
                light_range = light["light_range"]
            else:
                # Default value, nautical miles
                light_range = str(default_light_range)

            light_ini_file.write("Desc(" + str(i + 1) + ")=" + str(light) + "\n")

            # Define position, either through lat/long, or which buoy/mark it's associated with
            if light["light_buoy"] > 0:
                light_ini_file.write(
                    "Buoy(" + str(i + 1) + ")=" + str(light["light_buoy"]) + "\n"
                )
            else:
                light_ini_file.write(
                    "Long(" + str(i + 1) + ")=" + light["light_lon"] + "\n"
                )
                light_ini_file.write(
                    "Lat(" + str(i + 1) + ")=" + light["light_lat"] + "\n"
                )

            light_ini_file.write("Sequence(" + str(i + 1) + ")=" + sequence + "\n")
            light_ini_file.write("Range(" + str(i + 1) + ")=" + light_range + "\n")

            # For phase, use the node_index a bit like a hash input,
            # so lights on same node with the same length sequence are in phase
            phase_start = (light["node_index"] % len(sequence)) + 1
            light_ini_file.write(
                "PhaseStart(" + str(i + 1) + ")=" + str(phase_start) + "\n"
            )

            if light["light_height"] != "":
                # Note this is actually height above mean sea level from OSM definition,
                # might need adjusting to be above chart datum
                light_ini_file.write(
                    "Height(" + str(i + 1) + ")=" + light["light_height"] + "\n"
                )
                light_ini_file.write("Absolute(" + str(i + 1) + ")=1\n")
            else:
                light_ini_file.write(
                    "Height(" + str(i + 1) + ")=" + str(default_light_height) + "\n"
                )
                # Absolute=2 gives height relative to land, unless land height < 0, in which case height is relative to datum
                light_ini_file.write("Absolute(" + str(i + 1) + ")=2\n")

            if light["light_colour"] == "white":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "255" + "\n")
            elif light["light_colour"] == "red":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "0" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "0" + "\n")
            elif light["light_colour"] == "green":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "0" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "0" + "\n")
            elif light["light_colour"] == "blue":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "0" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "0" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "255" + "\n")
            elif light["light_colour"] == "yellow":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "0" + "\n")
            elif light["light_colour"] == "amber":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "128" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "0" + "\n")
            elif light["light_colour"] == "orange":  # Currently same as amber
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "128" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "0" + "\n")
            elif light["light_colour"] == "violet":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "128" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "255" + "\n")
            elif light["light_colour"] == "magenta":
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "128" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "0" + "\n")
            else:
                print("Light colour " + light["light_colour"] + " not mapped\n")
                light_ini_file.write("Red(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Green(" + str(i + 1) + ")=" + "255" + "\n")
                light_ini_file.write("Blue(" + str(i + 1) + ")=" + "255" + "\n")

            if (light["sector_start"] != "") and (light["sector_end"] != ""):
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

                light_ini_file.write(
                    "StartAngle(" + str(i + 1) + ")=" + str(start_angle) + "\n"
                )
                light_ini_file.write(
                    "EndAngle(" + str(i + 1) + ")=" + str(end_angle) + "\n"
                )
            else:
                light_ini_file.write("StartAngle(" + str(i + 1) + ")=" + "0" + "\n")
                light_ini_file.write("EndAngle(" + str(i + 1) + ")=" + "360" + "\n")

            light_ini_file.write("\n")
        else:
            print("Light character " + light["light_character"] + " not mapped\n")

    light_ini_file.close()

# Write a readme.txt file (For GMRT and OSM data)
readme_file_name = output_folder / "readme.txt"
readme_file = open(readme_file_name, "w")
readme_file.write(
    "Elevation and Bathymetry data from the Global Multi-Resolution Topography Synthesis (GMRT)\n"
)
readme_file.write(
    "For details, please see Ryan, W. B. F., S.M. Carbotte, J. Coplan, S. O'Hara, A. Melkonian, R. Arko,\n"
    "R.A. Weissel, V. Ferrini, A. Goodwillie, F. Nitsche, J. Bonczkowski, and R. Zemsky (2009),\n"
    "Global Multi-Resolution Topography (GMRT) synthesis data set, Geochem. Geophys. Geosyst., \n"
    "10, Q03014, doi:10.1029/2008GC002332.\n"
    "\n"
    "Coastline data and navigation aid data from from OpenStreetMap data, available\n"
    "under the Open Database License.\n"
    "\n"
    "NOT FOR USE IN NAVIGATION\n"
)
readme_file.close()

if use_osm2world:
    import osm_osmium_osm2world

    osm_osmium_osm2world.set_java_location(str(java_home.absolute()))

    print("Generating building models using OSM2World")

    model_path = output_folder / "Models"
    if not model_path.exists():
        model_path.mkdir()

    land_model_path = model_path / "LandObject"
    if not land_model_path.exists():
        land_model_path.mkdir()

    land_object_ini_file_name = output_folder / "landobject.ini"
    land_object_ini_file = open(land_object_ini_file_name, "w")

    land_object_ini_file.write(
        "Number = " + str(x_divisions_OSM2World * y_divisions_OSM2World) + "\n\n"
    )

    land_object_index = 0

    for x_division in range(x_divisions_OSM2World):
        print(str(int(100 * x_division / x_divisions_OSM2World)) + "% complete.")
        for y_division in range(y_divisions_OSM2World):

            land_object_index = land_object_index + 1

            min_lon = terrain_long + x_division * (
                terrain_long_extent / x_divisions_OSM2World
            )
            max_lon = terrain_long + (x_division + 1) * (
                terrain_long_extent / x_divisions_OSM2World
            )
            min_lat = terrain_lat + y_division * (
                terrain_lat_extent / y_divisions_OSM2World
            )
            max_lat = terrain_lat + (y_division + 1) * (
                terrain_lat_extent / y_divisions_OSM2World
            )

            model_name = "osm_" + str(x_division + 1) + "_" + str(y_division + 1)

            this_model_path = land_model_path / model_name
            if not this_model_path.exists():
                this_model_path.mkdir()

            osm_output_name = str(this_model_path / "filtered.osm")
            obj_output_name = str(this_model_path / "model.obj")

            # Clip from overall OSM file and make model
            osm_osmium_osm2world.clip_osm_file_osmium(
                osm_map_file, osm_output_name, min_lat, min_lon, max_lat, max_lon
            )

            # Check if the file size is small (less than 1kb), if so check if it only has 3 lines
            run_osm2world = True
            if os.path.getsize(osm_output_name) < 1024:
                with open(osm_output_name, "r") as osm_output_file:
                    osm_output_file_length = len(osm_output_file.readlines())
                if osm_output_file_length <= 3:
                    run_osm2world = False

            # Run OSM2World on the file
            if run_osm2world:
                osm_osmium_osm2world.osm_2_world(
                    str(osm_2_world_exe.absolute()), osm_output_name, obj_output_name
                )

            land_object_ini_file.write(
                "Type(" + str(land_object_index) + ") = " + model_name + "\n"
            )

            # Read model origin from the .obj file from a line like
            # '# Coordinate origin (0,0,0): lat 50.795517000000004, lon -1.1060862, ele 0'
            # Defaults in case we can't read the file:
            model_origin_lon = "0"
            model_origin_lat = "0"
            try:
                obj_file = open(obj_output_name, "r")
                for i in range(10):
                    # Assume the line we want is in the first 10 lines
                    obj_line = obj_file.readline()
                    if obj_line.startswith("# Coordinate origin"):
                        break
                obj_file.close()

                if obj_line != "":
                    obj_line_split = obj_line.split(" ")
                    model_origin_lat = obj_line_split[5].replace(",", "")
                    model_origin_lon = obj_line_split[7].replace(",", "")
            except FileNotFoundError:
                pass
                # Leave origin as default (0, 0)

            land_object_ini_file.write(
                "Lat(" + str(land_object_index) + ") = " + model_origin_lat + "\n"
            )
            land_object_ini_file.write(
                "Long(" + str(land_object_index) + ") = " + model_origin_lon + "\n"
            )
            land_object_ini_file.write(
                "HeightCorrection(" + str(land_object_index) + ") = " + "0.1" + "\n"
            )
            land_object_ini_file.write(
                "Rotation(" + str(land_object_index) + ") = " + "180" + "\n"
            )
            land_object_ini_file.write(
                "Morph(" + str(land_object_index) + ") = " + "1" + "\n"
            )

            model_ini_filename = this_model_path / "object.ini"
            model_ini_file = open(model_ini_filename, "w")
            model_ini_file.write('FileName="model.obj"\n')
            model_ini_file.write("ScaleFactor=1.0\n")
            model_ini_file.close()

    land_object_ini_file.close()

# End
print("Script end: " + str(datetime.datetime.now()))
