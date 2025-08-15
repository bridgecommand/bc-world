# Automating world model creation for Bridge Command 
This script is a work in progress, built to automate
the process of generating world models for Bridge 
Command from publicly available data.

This script uses [GMRT](https://www.gmrt.org/) 
bathymetry and elevation data, with 
[OpenStreetMap](https://www.openstreetmap.org/) 
data to increase detail of the coastline.

The OSM filtering sets a minimum height for 'land' 
areas, and a minimum depth for 'sea' areas, based
on the OpenStreetMap coastline data. If this is not
required, it can be disabled by changing the `use_osm`
variable.

To use the script, you will need to download the 
GMRT data for the area of interest in GeoTIFF format.
To do this, go to the [GMRT Map Tool](https://www.gmrt.org/GMRTMapTool/)
and select the area you want for your world model, 
then click on *Create Grid File*. Make sure the file
format is set to GeoTIFF, the *Mask* setting to 
*Unmasked*. You can choose the resolution, but will 
normally want *Maximum*. Download the file, and save
in the same location as the script.

You will also need to download the OpenStreetMap
coastlines shapefile from 
[https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip](https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip).
This is a large file (~800 Mb), and should be saved
in the same location as the script.

Make sure you have a Python environment available, 
and ensure you have the following packages installed:
- matplotlib
- numpy
- rasterio
- shapefile
- shapely

Open the script, and edit the variables in the top
section. The main ones are the path to the GMRT 
GeoTIFF file, and the output resolution.

**Warning**:
At the moment, the script is very slow to run, due
to the OSM filtering. To test the script, you may 
want to start with `use_osm = False`. As an example
, a model covering Portsmouth Harbour at 1025x1025 
resolution took about 5 hours to process on a modern
laptop.

When the script has finished, the files in the output
directory can be copied into a new folder in your
Bridge Command world models folder.

## Future work
It should be possible to increase the speed of the
OSM filtering significantly.

The terrain texture is currently very simplistic. 
This could be improved, either generating a more
realistic texture, or using available aerial or
satellite photography.

Additionally, from the OpenSeaMap project, there is
a lot of data on buoys and navigation lights and 
markers that could be used.

Ideally the script could also identify marina
features from OSM data, and generate a representation
of these automatically.