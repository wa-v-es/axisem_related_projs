import os
import yaml
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
import sys
matplotlib.rcParams.update({'font.size': 8})
from obspy import read
from obspy.core import Stream, Trace, UTCDateTime, Stats
from obspy.io.sac import SACTrace
###
# read event location
input_dir='/Users/keyser/Research/axisem/mars/s1222a_mars1D/input1D'
info_5deg = np.loadtxt(input_dir+'/mars.txt', dtype=str, skiprows=3)
st_dir = '/Users/keyser/Research/axisem/mars/seis_square10s'
####

with open(input_dir+'/inparam.source.yaml', 'r') as file:
    source_yaml = yaml.load(file, Loader=yaml.FullLoader)
loc_leaf = source_yaml['list_of_sources'][0]['s1222a_mars']['location']
event_latlon = loc_leaf['latitude_longitude']
event_depth = loc_leaf['depth']
#
# read time and displacement
station_key='11'
network_key='aux11'
time = np.loadtxt(st_dir + '/data_time.ascii')
disp1D = np.loadtxt(st_dir + '/%s.%s.ascii' % (network_key, station_key))


##Save to SAC after down-sampling

# create dir
os.makedirs(st_dir + '/sac', exist_ok=True)
# trace header
stats = Stats()
# stats.network='0'
stats.starttime = UTCDateTime(time[0])
stats.delta = UTCDateTime(time[1] - time[0])
stats.npts = len(time)

# sac header
sac_header = {}
sac_header['evla'] = event_latlon[0]
sac_header['evlo'] = event_latlon[1]
sac_header['evdp'] = float(event_depth)

# tr = Trace(data=disp1D[:, 2], header=stats)
# tr.plot()
# sys.exit()
# loop over stations
print('Saving to SAC...')
for ist, st in enumerate(info_5deg):
    print('%d / %d' % (ist + 1, len(info_5deg)), end='\r')
    # sac header
    sac_header['kstnm'] = st[0]
    sac_header['knetwk'] = st[1]
    sac_header['stla'] = float(st[2])
    sac_header['stlo'] = float(st[3])
    sac_header['stdp'] = float(st[5])
    # read data
    disp = np.loadtxt(st_dir + '/%s.%s.ascii' % (st[1], st[0]))
    # loop over channels
    for ich, ch in enumerate('RTZ'):
        # sac header
        sac_header['kcmpnm'] = ch
        # add sac header to trace header
        stats.sac = sac_header
        stats.station= st[0]
        stats.channel= ch
        # create and process trace
        tr = Trace(data=disp[:, ich], header=stats)
        tr.resample(20.)

        tr = tr.slice(UTCDateTime(0.), UTCDateTime(1800.))
        # create sac from trace
        sac = SACTrace.from_obspy_trace(tr)
        sac.write(st_dir + '/sac/%s.%s.%s.sac' % (st[1], st[0], ch))
    # sys.exit()
print('Done with %d stations.' % len(info_5deg))

sys.exit()
####################################
# draw a map of event and stations #
####################################
# plot


##
##Processing using obspy

plt.figure()
# draw map
# map = Basemap(projection='cyl', resolution='l', lon_0=180)
map = Basemap(height=16700000,width=12000000,
            resolution='l',area_thresh=700.,projection='omerc',\
            lon_0=-100,lat_0=15,lon_2=-170,lat_2=55,lon_1=-90,lat_1=-45)
map.drawcoastlines(linewidth=0.15)
map.fillcontinents(color='ivory',lake_color='lavenderblush')
map.drawmapboundary(fill_color='lavenderblush', linewidth=.1)
# draw event
map.scatter(event_latlon[1], event_latlon[0], latlon=True,
            s=150, c='salmon', marker='*', lw=0, zorder=100)
# draw stations
map.scatter(info_5deg[:, 3].astype(float), info_5deg[:, 2].astype(float), latlon=True,
            s=8, c='seagreen', edgecolor='black',linewidth=.1,marker='^', zorder=10,alpha=.55)
# plt.show()
plt.savefig('AK81.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
###
