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
input_dir='/Users/keyser/Research/axisem_related_projs/plumes/output_10sec_new_source/input'


plume_syn = '/Users/keyser/Research/axisem_related_projs/plumes/output_10sec_new_source/stations/100KM_sts/sac_bf_noise/'
no_plume_syn = '/Users/keyser/Research/axisem_related_projs/plumes/no_plume_10sec/output/stations/no_plume/sac_bf_noise/'


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
stations_subset = np.array([
    row for row in info_5deg
    if -3 <= float(row[2]) <= 3 and 6 <= float(row[3]) <= 15 and float(row[2])!=0
])

for row in stations_subset:
    row[0] = row[0].replace('la', 'L').replace('lo', 'L')
print('lenth of stations_subset=',len(stations_subset))
# print(stations_subset)
# sys.exit()
for ist, st in enumerate(stations_subset):
    print('%d / %d' % (ist + 1, len(info_5deg)), end='\r')
    st_name=st[0]
    # st_name='MCK'
    no_st_plume=read(no_plume_syn+'*{}.Z.sac'.format(st_name))
    st_plume=read(plume_syn+'*{}.*Z.sac'.format(st_name))
    no_st_plume[0].stats.station=st_name
    st_plume[0].stats.station=st_name

    # origin_time = get_sac_reftime(st_iris[0].stats.sac)
    # st_syn[0].stats.starttime=st_iris[0].stats.starttime


    #### to align the seismograms

    model = TauPyModel(model="ak135")
    arrivals = model.get_travel_times(source_depth_in_km=st_plume[0].stats.sac['evdp']/1000,distance_in_degree=st_plume[0].stats.sac['gcarc'],phase_list=["P"])
    arr_P=arrivals[0]

    # st_syn.trim(starttime=origin_time+arr_PP.time-400,endtime=origin_time+arr_PP.time+200)
    no_st_plume.filter('lowpass',freq=.04)
    st_plume.filter('lowpass',freq=.04)
    # st_iris[0].data=st_iris[0].data
    stream_both = Stream(traces=[no_st_plume[0], st_plume[0]])
    no_plume_all.append(no_st_plume[0])
    plume_all.append(st_plume[0])

    # stream_both.plot()
# sys.exit()
##
plt.ion()
###
#sort stream based on great circle distance
no_plume_all = Stream(sorted(no_plume_all, key=lambda tr: tr.stats.sac['gcarc']))
plume_all = Stream(sorted(plume_all, key=lambda tr: tr.stats.sac['gcarc']))

#
fig, ax = plt.subplots(1,3, dpi=150,figsize=(15.5, 7))
N_old=0
origin_time = get_sac_reftime(no_plume_all[0].stats.sac)
for j in range(3):
    # print(j)
    N_new = (j + 1) * 17

    if N_new >= len(no_plume_all):  # last bit
        st_temp = no_plume_all[N_old:]
    else:
        st_temp = no_plume_all[N_old:N_new]  # current chunk
    N_old = N_new
    print(st_temp,'\n')
    print('################################')

    ax[j].set_ylim(18, 0)
    # ax[j].set_xlim(650,1195)
    ax[j].set_yticks(np.linspace(0,18,19))
    ax[j].tick_params(axis='y',left=False,pad=1)
    # ax[j].yaxis.labelpad = -12
    ax[j].grid(which='major', axis='x',color='dimGrey', linestyle='--',linewidth=.5,alpha=.95)
    ax[j].xaxis.set_minor_locator(MultipleLocator(50))
    ax[j].xaxis.set_major_locator(MultipleLocator(100))
    # ax[j].annotate('Vertical ACF AEB09', xy=(.6, 1.05), xycoords='axes fraction')
    # stack_ax1=iris_all[:17]
    for i in range(len(st_temp)):
        auto=st_temp[i]
        model = TauPyModel(model="ak135")
        # for synthetic, dep is in km. In iris, it is in mt.
        arrivals = model.get_travel_times(source_depth_in_km=auto.stats.sac['evdp']/1000,distance_in_degree=auto.stats.sac['gcarc'],phase_list=("P","pP","sP",'pS','S'))
        arr_P=arrivals[0]
        # try:
        #     arrivals = model.get_travel_times(source_depth_in_km=auto.stats.sac['evdp']/1000,distance_in_degree=auto.stats.sac['gcarc'],phase_list=["P"])
        #     # arrivals = model.get_travel_times(source_depth_in_km=st_iris[0].stats.sac['evdp']/1000,distance_in_degree=st_iris[0].stats.sac['gcarc'],phase_list=["P"])
        #     arr_P=arrivals[0]
        #     print(arr_P.time,'sec\n')
        # except:
        #     print("P shadow zone reached!\n")

        auto.trim(starttime=origin_time+arr_P.time-20)#,endtime=origin_time+arr_P.time+160)

        if auto.stats.npts == 0:
            print('No samples in trace\n')
            continue

        time = np.arange(auto.stats.npts) * auto.stats.delta
        # time=time+arr_PP.time-400 # shifts onset to 0 sec
        auto.data /= np.max(np.abs(auto.data))
        l1,=ax[j].plot(time,auto.data + i+1,  lw=0.5, color='teal',label='real')
        for arr in arrivals:
            ax[j].scatter(arr.time - arr_P.time+20, i+1,s=10,marker='o',facecolor='navy', edgecolor='black',alpha=.95,linewidth=.15,zorder=10)

        ####
        ##plotting synthetic on top
        auto_s=plume_all.select(station=auto.stats.station)[0]
        auto_s.trim(starttime=origin_time+arr_P.time-20)#,endtime=origin_time+arr_P.time+160)

        time_s = np.arange(auto_s.stats.npts) * auto_s.stats.delta
        ### ADDED 10 sec in synthetic!!!!!!
        # time_s=time_s+arr_PP.time-400+10 # shifts onset to 0 sec
        # time_s=time_s+arr_PP.time-400

        auto_s.data /= np.max(np.abs(auto_s.data))

        l2,=ax[j].plot(time_s,auto_s.data + i+1,  lw=0.55,ls='--', color='brown',label='Axisem')


        # ax[j].fill_between(time, i+1, auto.data + i+1, lw=0.55,
        #                   color=cm.viridis(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
        # plt.setp(auto.stats.station,fontsize=3)
    #
    st_label=[' ']
    for i in range(len(st_temp)):
        auto=st_temp[i]
        reso_st=[]
        # reso_filter.append('{:.2f}, {:.2f}, {:.2f}'.format(auto.stats.tpdelay,auto.stats.Dt_mean,auto.stats.Dt_std))
        st_label.append(auto.stats.station)
    plt.setp(ax[j].get_yticklabels(),fontsize=7)
    ax[j].yaxis.set_major_formatter(ticker.FixedFormatter(st_label))
# ax[1].set_xlabel('Time after eq (s)')
plt.legend([l1,l2],['No plume','Plume'], loc=[-1,1.02],ncol=3,fontsize=9,handletextpad=.5,borderaxespad=1.5,columnspacing=1)
fig.text(0.52, 0.03, 'Centered around P-20 lowpass (0.04 Hz/25 s)',fontsize=11, ha='center', va='center')
# plt.show()
#
# plt.savefig('noise_plume_noPlume_10mesh_Z_25s.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
sys.exit()
############

model = TauPyModel(model="ak135")
arrivals = model.get_travel_times(source_depth_in_km=st_iris[0].stats.sac['evdp']/1000,distance_in_degree=st_iris[0].stats.sac['gcarc'],phase_list=['P','pP','sP',"PP"])
plt.ion()

# plot synthetic with P, S, and other phases.
fig, ax = plt.subplots(2, sharex=True, dpi=150)
# for ich, ch in enumerate('RTZ'):
#     # change unit to mm
ax[j].xaxis_date()
ax[1].xaxis_date()
ax[j].plot(st_iris[0].times("utcdatetime"), st_iris[0].data, lw=.8,c='cadetblue', label='real')
ax[1].plot(st_syn[0].times("utcdatetime"), st_syn[0].data, lw=.8,c='indianred' ,label='syn')
for i in range(len(arrivals)):
    ax[1].scatter(origin_time+arrivals[i].time,0,s=17,marker='d',facecolor='rebeccapurple', edgecolor='black',alpha=.5,linewidth=.75,zorder=10)
    ax[j].scatter(origin_time+arrivals[i].time,0,s=17,marker='d',facecolor='rebeccapurple', edgecolor='black',alpha=.5,linewidth=.75,zorder=10)

#

fig.autofmt_xdate()
plt.legend()
plt.savefig('AK81_real_syn.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()

# ax[1].text(.95, .2, 'channel = ' + ch, transform = ax[ich].transAxes, ha='right', va='top')
# ax[1].set_ylabel('Amplitude (mm)')
# ax[j].set_xlim(time[0], time[-1])
# plt.xlabel('Time after source origin (s)')
##
#cd '/Users/keyser/Research/axisem/forPC_output/stations/5deg_array/sac'
st=read('*long80.Z.sac')
arrivals = model.get_travel_times(source_depth_in_km=600,distance_in_degree=80,phase_list=['P','pP','sP',"PP"])
fig = plt.figure(figsize=(11,8))
st.plot(color='cadetblue',bgcolor='ivory', fig=fig, orientation="horizontal")
for i in range(len(arrivals)):
    plt.scatter(st[0].stats.starttime+arrivals[i].time,0,s=37,marker='o',facecolor='rebeccapurple',alpha=.95,zorder=10)
plt.show()
