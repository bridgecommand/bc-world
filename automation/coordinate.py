import math

def deg2num(lat_deg, lon_deg, zoom):
  lat_rad = math.radians(lat_deg)
  n = 2.0 ** zoom
  xtile = int((lon_deg + 180.0) / 360.0 * n)
  ytile = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
  return (zoom,xtile, ytile)

def num2deg(xtile, ytile, zoom):
  n = 2.0 ** zoom
  lon_deg = xtile / n * 360.0 - 180.0
  lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
  lat_deg = math.degrees(lat_rad)
  return (lat_deg, lon_deg)

# Read API key from apiKey.txt
f = open("apiKey.txt", "r")
apiKey = f.readline()


lowerLeftLongRequest = 50.7965
lowerLeftLatRequest = -1.1115
zoomRequest = 14

# Find tile containing the requested point
(zoom, xtile, ytile) = deg2num(lowerLeftLongRequest,lowerLeftLatRequest,zoomRequest)
tileURL = "https://tile.nextzen.org/tilezen/terrain/v1/512/terrarium/"+str(zoom)+"/"+str(xtile)+"/"+str(ytile)+".png?api_key="+str(apikey)

# Find lat/long information about the corners. We will get the NW point for each tile. Adding 1 to tile number gives next tile in the SE direction
(maxLat, minLong) = num2deg(xtile, ytile, zoom)
(minLat, maxLong) = num2deg(xtile+1, ytile+1, zoom)

print("Zoom.      : ", zoom)
print("X tile     : ", xtile)
print("Y tile     : ", ytile)
print("Long.      : ", minLong)
print("Lat.       : ", minLat)
print("Long extent: ", maxLong - minLong)
print("Lat extent : ", maxLat  - minLat)
print(tileURL)