import gpxpy
import numpy as np
import scipy.io
import math
from geopy.distance import vincenty
from geopy.distance import great_circle
import datetime
import matplotlib.pyplot as plt

filename = '../GPX/leif_324.gpx'
filename_out = 'leif_324.mat'
gps = gpxpy.parse(open(filename)) #'activity_1679903652.gpx'))
points = gps.tracks[0].segments[0].points
segment = gps.tracks[0].segments[0]



lenRun = len(points)
print('points in run:', lenRun)

lat = np.empty(0)
lon = np.empty(0)
alt = np.empty(0)


for point in points:
	lat = np.append(lat,point.latitude)
	lon = np.append(lon,point.longitude)
	alt = np.append(alt,point.elevation)


try:
	test = points[0].time
	time_run = np.empty(0)
	for point in points:
		time_run = np.append(time_run,point.time)
	#time to sec
	Ts = np.array([])
	for i, stamp in enumerate(time_run[0:-1]):
		dt = time_run[i+1]-time_run[i]
		dt = dt.seconds
		Ts = np.append(Ts,dt)
	Ts = Ts[::-1]
except TypeError:
	print('No Time Information')



# # reverse gpx data to start at NW Thurman
# lat = lat[::-1]
# lon = lon[::-1]
# alt = alt[::-1]


lon_m = np.empty(0)
lat_m = np.empty(0)
alt_m = alt*0.3048


for n in range(0,lenRun):
	# longitude
	l0 = (lon[n],0)
	l1 = (0,0)
	crow = vincenty(l0,l1).meters
	lon_m = np.append(lon_m,crow)
	# latitude
	l0 = (lat[n],0)
	l1 = (0,0)
	crow = vincenty(l0,l1).meters
	lat_m = np.append(lat_m,crow)

plt.plot(lon_m,lat_m)
plt.show()


# 	d0 = (lon[n],lat[n])
# 	d1 = (lon[n+1],lat[n+1])
# 	crow_fly = vincenty(d0, d1).meters
# 	# crow_fly = great_circle(d0, d1).meters
# 	wolf_run = math.sqrt(crow_fly**2 + alt_run_d[n]**2)
# 	dist_wolf_run = np.append(dist_wolf_run,wolf_run)
# 	dist_crow_fly = np.append(dist_crow_fly,crow_fly)
# 	l0 = (lon[n],0)
# 	l1 = (lon[n+1],0)
# 	crow = vincenty(l0,l1).meters
# 	wolf = math.sqrt(crow**2 + alt_run_d[n]**2)
# 	lon_dist_wolf_run = np.append(lon_dist_wolf_run,wolf)
# 	l0 = (0,lat[n])
# 	l1 = (0,lat[n+1])
# 	crow = vincenty(l0,l1).meters
# 	wolf = math.sqrt(crow**2 + alt_run_d[n]**2)
# 	lat_dist_wolf_run = np.append(lat_dist_wolf_run,wolf)

# variance = np.var(wolf_run)
# print(variance)

scipy.io.savemat(filename_out,mdict={'lon_m': lon_m, 'lat_m':lat_m,'alt_m':alt_m,'time':Ts})



