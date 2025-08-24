## this script compares the synthetic and real seismograms.
# synthetic computed using axisem3d
###
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
from obspy.io.sac.util import get_sac_reftime
from obspy.taup import TauPyModel
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
########
input_dir='/Users/keyser/Research/axisem_related_projs/plumes/plumes_iaspi91_10sec_new_loc_with_wave/input'


plume_syn = '/Users/keyser/Research/axisem_related_projs/plumes/plumes_iaspi91_10sec_new_loc_with_wave/simu1D/output/stations/plume_2lat/sac_bf/'
no_plume_syn = '/Users/keyser/Research/axisem_related_projs/plumes/plumes_iaspi91_10sec_new_loc_wave_disabled/simu1D/output/stations/plume_2lat/sac_bf/'

info_5deg = np.loadtxt(input_dir+'/grid_stations.txt', dtype=str, skiprows=3)

with open(input_dir+'/inparam.source.yaml', 'r') as file:
    source_yaml = yaml.load(file, Loader=yaml.FullLoader)
loc_leaf = source_yaml['list_of_sources'][0]['the_only_source']['location']
event_latlon = loc_leaf['latitude_longitude']
event_depth = loc_leaf['depth']
###
# starttime=ori.time+arr_PP.time-400 for st_iris
plume_all=Stream()
no_plume_all=Stream()

# row 2 is lat. row 3 is long.
lats=(5,15)
# lats=(5,)
lons=(5,15)
# lons=(5,)
stations_subset = np.array([
    row for row in info_5deg
    if float(row[2]) in lats and float(row[3]) in lons
])
# if float(row[2]) == 15 and float(row[3]) == 15
for row in stations_subset:
    row[0] = row[0].replace('la', '').replace('lo', '')
    # row[0] = row[0].replace('L', '')

stations_subset=np.delete(stations_subset, 2,axis=0)
# to plot just one trace
# stations_subset=np.delete(stations_subset, 1,axis=0)
# stations_subset=np.delete(stations_subset, 1,axis=0)

print('lenth of stations_subset=',len(stations_subset))
# print(stations_subset)
# sys.exit()
for ist, st in enumerate(stations_subset):
    print('%d / %d' % (ist + 1, len(info_5deg)), end='\r')
    st_name=st[0]
    # st_name='MCK'
    no_st_plume=read(no_plume_syn+'{}.R.sac'.format(st_name))
    st_plume=read(plume_syn+'{}.R.sac'.format(st_name))
    no_st_plume[0].stats.station=st_name
    st_plume[0].stats.station=st_name

    # origin_time = get_sac_reftime(st_iris[0].stats.sac)
    # st_syn[0].stats.starttime=st_iris[0].stats.starttime


    #### to align the seismograms

    model = TauPyModel(model="ak135")
    arrivals = model.get_travel_times(source_depth_in_km=st_plume[0].stats.sac['evdp']/1000,distance_in_degree=st_plume[0].stats.sac['gcarc'],phase_list=["P"])
    arr_P=arrivals[0]

    # st_syn.trim(starttime=origin_time+arr_PP.time-400,endtime=origin_time+arr_PP.time+200)
    no_st_plume.filter('lowpass',freq=.2)
    st_plume.filter('lowpass',freq=.2)
    # st_iris[0].data=st_iris[0].data
    stream_both = Stream(traces=[no_st_plume[0], st_plume[0]])
    no_plume_all.append(no_st_plume[0])
    plume_all.append(st_plume[0])

    # stream_both.plot()
# sys.exit()
##
# plt.ion()
###
#sort stream based on great circle distance
no_plume_all = Stream(sorted(no_plume_all, key=lambda tr: tr.stats.sac['gcarc']))
plume_all = Stream(sorted(plume_all, key=lambda tr: tr.stats.sac['gcarc']))

#
plt.ion()
plt.rcParams['font.size'] = 12

fig, ax = plt.subplots(1,1, dpi=150,figsize=(12, 4))

origin_time = get_sac_reftime(no_plume_all[0].stats.sac)


st_temp = no_plume_all  # current chunk
# N_old = N_new
print(st_temp,'\n')
print('################################')

ax.set_ylim(4, 0)
# ax.set_ylim(2, 0)

ax.set_xlim(0,450)
ax.set_yticks(np.linspace(0,4,5))
ax.tick_params(axis='y',left=False,pad=1)
# ax.yaxis.labelpad = -12
ax.grid(which='major', axis='x',color='dimGrey', linestyle='--',linewidth=.5,alpha=.95)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.xaxis.set_major_locator(MultipleLocator(50))
# ax.annotate('Vertical ACF AEB09', xy=(.6, 1.05), xycoords='axes fraction')
# stack_ax1=iris_all[:17]
for i,auto in enumerate(st_temp):
    auto=st_temp[i]
    model = TauPyModel(model="ak135")
    # for synthetic, dep is in km. In iris, it is in mt.
    arrivals = model.get_travel_times(source_depth_in_km=auto.stats.sac['evdp']/1000,distance_in_degree=auto.stats.sac['gcarc'],phase_list=("P","pP","sP",'pS','S','SP','PS'))
    arr_P=arrivals[0]
    # try:
    #     arrivals = model.get_travel_times(source_depth_in_km=auto.stats.sac['evdp']/1000,distance_in_degree=auto.stats.sac['gcarc'],phase_list=["P"])
    #     # arrivals = model.get_travel_times(source_depth_in_km=st_iris[0].stats.sac['evdp']/1000,distance_in_degree=st_iris[0].stats.sac['gcarc'],phase_list=["P"])
    #     arr_P=arrivals[0]
    #     print(arr_P.time,'sec\n')
    # except:
    #     print("P shadow zone reached!\n")

    auto.trim(starttime=origin_time+arr_P.time-50)#,endtime=origin_time+arr_P.time+300)

    if auto.stats.npts == 0:
        print('No samples in trace\n')
        continue

    time = np.arange(auto.stats.npts) * auto.stats.delta
    # time=time+arr_PP.time-400 # shifts onset to 0 sec
    auto.data /= np.max(np.abs(auto.data))
    l1,=ax.plot(time,auto.data*.6 + i+1,  lw=0.95, color='teal',ls='--',label='real')
    for arr in arrivals:
        ax.scatter(arr.time - arr_P.time+50, i+1,s=10,marker='o',facecolor='navy', edgecolor='black',alpha=.95,linewidth=.15,zorder=10)
        ax.text(arr.time - arr_P.time+50, i+.5, arr.name, bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1.5},fontsize=9.5,c='navy', rotation='vertical',ha='center')


    ####
    ##plotting synthetic on top
    auto_s=plume_all.select(station=auto.stats.station)[0]
    auto_s.trim(starttime=origin_time+arr_P.time-50)#,endtime=origin_time+arr_P.time+300)

    time_s = np.arange(auto_s.stats.npts) * auto_s.stats.delta
    ### ADDED 10 sec in synthetic!!!!!!
    # time_s=time_s+arr_PP.time-400+10 # shifts onset to 0 sec
    # time_s=time_s+arr_PP.time-400

    auto_s.data /= np.max(np.abs(auto_s.data))

    l2,=ax.plot(time_s,auto_s.data*.6 + i+1,  lw=0.85, color='brown',label='Axisem')


    # ax.fill_between(time, i+1, auto.data + i+1, lw=0.55,
    #                   color=cm.viridis(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
    # plt.setp(auto.stats.station,fontsize=3)
#
st_label=[' ']
for i in range(len(st_temp)):
    auto=st_temp[i]
    reso_st=[]
    # reso_filter.append('{:.2f}, {:.2f}, {:.2f}'.format(auto.stats.tpdelay,auto.stats.Dt_mean,auto.stats.Dt_std))
    st_label.append(auto.stats.station)
plt.setp(ax.get_yticklabels(),fontsize=7)
ax.yaxis.set_major_formatter(ticker.FixedFormatter(st_label))

# ax[1].set_xlabel('Time after eq (s)')
plt.legend([l1,l2],['No plume','Plume'], loc=[.65,1.02],ncol=3,fontsize=12,handletextpad=.5,borderaxespad=1.5,columnspacing=1)
fig.text(0.42, .92, 'Centered around P-20 lowpass (0.1 Hz )',fontsize=12, ha='center', va='center')
fig.text(0.5, .01, 'Time (s)',fontsize=12, ha='center', va='center')

plt.show()
#
# plt.savefig('3seis_noise_plume_noPlume_10mesh_Z_10s.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
sys.exit()
############
