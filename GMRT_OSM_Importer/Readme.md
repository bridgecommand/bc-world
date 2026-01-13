# Automating world model creation for Bridge Command 
This script is a work in progress, built to automate
the process of generating world models for Bridge 
Command from publicly available data.

This script uses [GMRT](https://www.gmrt.org/) 
bathymetry and elevation data, with 
[OpenStreetMap](https://www.openstreetmap.org/) 
data to increase detail of the coastline.

Optionally OSM2World is used to generate models of
buildings and land features. This requires Java, 
OSM2World and Osmosis to be installed, and requires  
the alpha version of Bridge Command (5.10.4-alpha.2 
onwards) which supports 'morphing' of models to fit 
the terrain shape.

The OSM filtering sets a minimum height for 'land' 
areas, and a minimum depth for 'sea' areas, based
on the OpenStreetMap coastline data. If this is not
required, it can be disabled by changing the 
`use_osm_coastline` variable.

The script also uses additional data from OpenStreetMap 
to add buoys and other aids to navigation. This is 
currently in an early stage of development, and imports 
lateral and cardinal buoys, lateral and special posts,
and mooring buoys. Importing navigation lights is currently 
in the early stages of development. The OpenStreetMap for the relevant
area is downloaded automatically. If importing this
data is not required, it can be disabled by changing
the `use_osm_map` variable.

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
- pyshp (shapefile)
- shapely

And if you want to generate building models automatically:
- osmium
- You will also need a Java JRE or JDK available, and set the java_home variable in main.py to point to this
- You will need to download OSM2World from https://osm2world.org/download/files/latest/OSM2World-latest-bin.zip (or the 0.5 version once this has been released) and set the osm_2_world_exe variable to point to this

Open the script, and edit the variables in the top
section. The main ones are the path to the GMRT 
GeoTIFF file, and the output resolution.

If you want to use the use_osm2world option to enable
the generation of building models, a java runtime must
be available, you will need to edit the variable java_home
to point to this.

When the script has finished, the files in the output
directory can be copied into a new folder in your
Bridge Command world models folder.

## Future work
The terrain texture is currently very simplistic. 
This could be improved, either generating a more
realistic texture, or using available aerial or
satellite photography.

Only a subset of buoy and fixed navigation mark
information, and a subset of navigation light information
is currently loaded.

Ideally the script could also identify marina
features from OSM data, and generate a representation
of these automatically.