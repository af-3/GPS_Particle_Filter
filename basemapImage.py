from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import gpxpy
import numpy as np
from scipy.interpolate import spline
import scipy.io


gpx_leif = gpxpy.parse(open('/Users/afockele2/Google Drive/School/PSU/2017_3_Summer/Project/GPX/nw-leif-erikson-drive-trail.gpx'))
points_leif = gpx_leif.tracks[0].segments[0].points
lat_leif = np.empty(0)
lon_leif = np.empty(0)
alt_leif = np.empty(0)

for point in points_leif:
	lat_leif = np.append(lat_leif,point.latitude)
	lon_leif = np.append(lon_leif,point.longitude)
	alt_leif = np.append(alt_leif,point.elevation)
# reverse gpx data to start at NW Thurman
lat_leif = lat_leif[::-1]
lon_leif = lon_leif[::-1]
alt_leif = alt_leif[::-1]

# scipy.io.savemat('leif_cords.mat',mdict={'lon':lon_leif,'lat':lat_leif,'alt':alt_leif})

# OVERLAY BASEMAP
map = plt.figure(2)
map = Basemap(llcrnrlon=-122.80,llcrnrlat=45.53,urcrnrlon=-122.72,urcrnrlat=45.60, epsg=4326, projection = 'tmerc')
#http://server.arcgisonline.com/arcgis/rest/services
# map.shadedrelief()
map.plot(lon_leif,lat_leif,color = 'green',lw = 1.0)
map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
plt.savefig('../Reports/basemap.png')
plt.show()
