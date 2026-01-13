import subprocess
import os
import osmium

def set_java_location(java_home):
    os.environ["JAVA_HOME"]=java_home


def clip_osm_file_osmium(input_file, output_file, min_lat, min_lon, max_lat, max_lon):
    # Iterate first to find relations to add (needed to contain building parts)
    fp = osmium.FileProcessor(input_file).with_locations().with_areas().with_filter(
        osmium.filter.KeyFilter('building', 'building:part'))
    relations_to_add = []
    for obj in fp:
        if obj.is_area():
            lon_list = []
            lat_list = []
            for outer in obj.outer_rings():
                for n in outer:
                    if n.location.valid():
                        lon_list.append(n.lon)
                        lat_list.append(n.lat)
            if len(lat_list) > 0 and len(lon_list) > 0:
                mean_lat = sum(lat_list) / len(lat_list)
                mean_lon = sum(lon_list) / len(lon_list)
                if min_lat < mean_lat <= max_lat:
                    if min_lon < mean_lon <= max_lon:
                        # Only add if the simple average of the nodal positions is within the required bounding box
                        relations_to_add.append(obj.orig_id())

    # Main iteration to add ways, plus the relations we found above
    fp = osmium.FileProcessor(input_file).with_locations().with_filter(
        osmium.filter.KeyFilter('building', 'building:part'))
    with (osmium.BackReferenceWriter(output_file, ref_src=input_file,
                                     overwrite=True) as writer):
        for obj in fp:
            if obj.is_way():
                lon_list = []
                lat_list = []
                for n in obj.nodes:
                    if n.location.valid():
                        lon_list.append(n.lon)
                        lat_list.append(n.lat)
                if len(lat_list) > 0 and len(lon_list) > 0:
                    mean_lat = sum(lat_list) / len(lat_list)
                    mean_lon = sum(lon_list) / len(lon_list)
                    if min_lat < mean_lat <= max_lat:
                        if min_lon < mean_lon <= max_lon:
                            # Only add if the simple average of the nodal positions is within the required bounding box
                            writer.add(obj)
            if obj.is_relation():
                if obj.id in relations_to_add:
                    writer.add(obj)

    add_bounding_box = False

    if add_bounding_box:
        # Manually add bounding box line, like '<bounds minlat="50.7585000" minlon="-1.1780000" maxlat="50.8519000" maxlon="-1.0583000"/>'
        # as third line
        f = open(output_file, 'r', encoding='UTF-8')
        raw_lines = f.readlines()
        f.close()

        if len(raw_lines) >= 3:
            f = open(output_file, 'w', encoding='UTF-8')
            f.write(raw_lines[0])
            f.write(raw_lines[1])
            bounding_string = f'<bounds minlat="{min_lat}" minlon="{min_lon}" maxlat="{max_lat}" maxlon="{max_lon}"/>\n'
            f.write(bounding_string)
            for i in range(2, len(raw_lines)):
                f.write(raw_lines[i])
            f.close()


def osm_2_world(path_to_osm2world, input_file, output_file):
    subprocess.run(args=[path_to_osm2world,
                         "convert",
                         "--lod=4",
                         "--input=" + input_file,
                         "--output=" + output_file])
