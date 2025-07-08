import os
import yaml
import numpy as np
from obspy.taup import TauPyModel
from obspy.taup.taup_geo import calc_dist,calc_dist_azi
from obspy.geodetics import gps2dist_azimuth
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib
import sys
matplotlib.rcParams.update({'font.size': 8})
from obspy import read
from obspy.core import Stream, Trace, UTCDateTime, Stats
from obspy.io.sac import SACTrace
#####
###
# read event location
input_dir='/Users/keyser/Research/axisem_related_projs/plumes/input'


info_arr = np.loadtxt(input_dir+'/stations_new.txt', dtype=str, skiprows=3)

# st_dir = '/Users/keyser/Research/axisem/loyalty_isl/output_10sec_2HD/stations/AK_81'
st_dir = '/Users/keyser/Research/axisem_related_projs/plumes/output/stations/100KM_sts'

####

with open(input_dir+'/inparam.source.yaml', 'r') as file:
    source_yaml = yaml.load(file, Loader=yaml.FullLoader)
loc_leaf = source_yaml['list_of_sources'][0]['the_only_source']['location']
event_latlon = loc_leaf['latitude_longitude']
event_depth = loc_leaf['depth']
#
# read time and displacement
# station_key='BAE'
network_key='II'
time = np.loadtxt(st_dir + '/data_time.ascii')
# disp1D = np.loadtxt(st_dir + '/%s.%s.ascii' % (network_key, station_key))

# sys.exit()
##Save to SAC after down-sampling

# create dir
os.makedirs(st_dir + '/sac_bf', exist_ok=True)
# os.makedirs(st_dir + '/sac_bf_noise', exist_ok=True)

# trace header
stats = Stats()
stats.network='II'
stats.starttime = UTCDateTime(time[0])
stats.delta = UTCDateTime(time[1] - time[0])
stats.npts = len(time)

# sac header
sac_header = {}
sac_header['evla'] = event_latlon[0]
sac_header['evlo'] = event_latlon[1]
sac_header['evdp'] = float(event_depth)
###
# sys.exit()
# loop over stations
add_noise=False
print('Saving to SAC...')
for ist, st in enumerate(info_arr):
    print('%d / %d' % (ist + 1, len(info_arr)), end='\r')
    # get PP arrival
    # model = TauPyModel(model="ak135")
    # dist=calc_dist(event_latlon[0],event_latlon[1],float(st[2]),float(st[3]),6400,0)
    dist_azi_bazi=calc_dist_azi(event_latlon[0],event_latlon[1],float(st[2]),float(st[3]),6400,0)
    #
    # arrivals = model.get_travel_times(source_depth_in_km=float(event_depth)/1000,distance_in_degree=dist,phase_list=["PP"])
    # arr_PP=arrivals[0]
    # try:
    #     arr_pP=model.get_travel_times(source_depth_in_km=float(event_depth)/1000,distance_in_degree=dist,phase_list=["pP"])[0]
    # except:
    #     print('no pP :/')
    # arris_P = model.get_travel_times(source_depth_in_km=float(event_depth)/1000,distance_in_degree=dist,phase_list=["P"])
    # try:
    #     arr_P=arris_P[0]
    # except:
    #     print('no P')
    # starttime=stats.starttime+arr_PP.time-400
    # endtime=stats.starttime+arr_PP.time+200
    # sys.exit()
    # sac header
    sac_header['b'] = 0
    sac_header['kstnm'] = st[0]
    sac_header['knetwk'] = st[1]
    sac_header['stla'] = float(st[2])
    sac_header['stlo'] = float(st[3])
    sac_header['stdp'] = float(st[5])
    sac_header['stel'] = float(st[4])
    sac_header['gcarc'] = dist_azi_bazi[0]
    sac_header['dist'] = dist_azi_bazi[0]*111
    sac_header['az'] = dist_azi_bazi[1]
    sac_header['baz'] = dist_azi_bazi[2]

    # read data
    disp = np.loadtxt(st_dir + '/%s.%s.ascii' % (st[1], st[0]))
    # sys.exit()
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
        tr.filter('bandpass',freqmin=.05,freqmax=.5)
        # tr.resample(20.)
        # tr.trim(starttime,endtime)
        # tr = tr.slice(UTCDateTime(0.), UTCDateTime(1800.))
        # create sac from trace
        sac = SACTrace.from_obspy_trace(tr)
        sac.write(st_dir + '/sac_bf/%s.%s.%s.sac' % (st[1], st[0], ch))

        if add_noise:
            max_amplitude = np.max(np.abs(tr.data))
            noise_std = 0.05 * max_amplitude

            # Gaussian noise
            noise = np.random.normal(0, noise_std, len(tr.data))
            # Add the noise to the trace
            tr.data += noise
            sac = SACTrace.from_obspy_trace(tr)
            sac.write(st_dir + '/sac_bf_noise/%s.%s.%s.sac' % (st[1], st[0], ch))

    # sys.exit()
print('Done with %d stations.' % len(info_arr))

sys.exit()
# copy binary
#cp axisem3d ./simu1D/
####################################
# check arrivals wrt synthetics
####################################
# plot

##
model = TauPyModel(model="ak135")
# dist_azi_bazi=calc_dist_azi(-21.19,170.27,61.1416, -148.1751,6400,0)
print(dist_azi_bazi)
arrivals = model.get_travel_times(source_depth_in_km=float(event_depth)/1000,distance_in_degree=dist,phase_list=['P','pP','sP',"PP"])
plt.ion()

# plot synthetic with P, S, and other phases.
fig, ax = plt.subplots(1, sharex=True, dpi=150)
# for ich, ch in enumerate('RTZ'):
#     # change unit to mm
ax.xaxis_date()
# ax[1].xaxis_date()
ax.plot(tr.times("utcdatetime"), tr.data, lw=.8,c='cadetblue', label='real')
# ax[1].plot(st_syn[0].times("utcdatetime"), st_syn[0].data, lw=.8,c='indianred' ,label='syn')
for i in range(len(arrivals)):
    # ax[1].scatter(origin_time+arrivals[i].time,0,s=17,marker='d',facecolor='rebeccapurple', edgecolor='black',alpha=.5,linewidth=.75,zorder=10)
    ax.scatter(stats.starttime+arrivals[i].time,0,s=17,marker='d',facecolor='rebeccapurple', edgecolor='black',alpha=.5,linewidth=.75,zorder=10)

#

fig.autofmt_xdate()
plt.legend()
plt.show()
