# reads moho_raw_datafile and comverts it to a smooth .nc file. Also, plots the point data and .nc as meshes.
# from https://github.com/benjaminfernando/axisem3d-3dmodels/tree/main/template_input/python%20scripts.
# original is axisem_moho_ori.py

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter, zoom

import miller_alaskamoho_srl2018 as alaskamoho
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os.path as path
import stripy as stripy
import sys
from scipy.spatial import cKDTree
import cmcrameri as cm
#

# moho_raw_datafile = 'depthtomoho.xyz'

moho_raw_datafile = 'MM_c1_tr2.XYZ'
moho_raw_datafile = 'dipping_moho/dip_prll_N_c1_tr2.XYZ'



##
# mynorm = plt.Normalize(vmin=-20000, vmax=20000)
mynorm = plt.Normalize(vmin=-30000, vmax=30000)

# sm = plt.cm.ScalarMappable(cmap='cmc.vik', norm=mynorm)


def load_moho_data(plot=False):
    data = np.hsplit(np.loadtxt(moho_raw_datafile),3)

    latitude = np.linspace(-89.999,89.999,num=180) #don't flip this coordinate as else not ascendingly sorted
    longitude = np.linspace(-179.999,179.999,num=360)

    depth = np.flip(np.reshape(data[2],(180,360)),0)

    if plot ==True:
        plot_moho(latitude, longitude, depth)

    return latitude, longitude, depth

def write_data_to_netcdf_file(smooth_overall=True, smooth_poles=False,resample=True, sigma=1, resample_factor =10, write_file = 'moho_rf.nc'):
    latitude_points, longitude_points, depth_points = load_moho_data()
    depth_points = depth_points*1000 #convert to metres here

    ###I THOUGHT SMOOTHING AT POLES MIGHT BE NEEDED; ROUTINE BELOW FOR LINEAR INTERPOLATION BETWEEN 80/90 N/S - DOESN'T FIX ISSUES BUT HERE
    ### IN CASE SOMEONE ELSE WANTS TO TRY IT. CAN ALSO SET ALL POINTS AT 90 DEGREES TO BE THE SAME: THIS IS AN ARTEFACT OF GOING FROM GRID-CENTRED
    ###TO CELL CENTRED (POINTS AT 89.5 DEGREES MOVE TO 90 DEGREES).

    #set points at poles to smooth?
    mean_north_pole_moho_depth = np.mean(depth_points[-1,:]) #origin is in the South-West corner as coordinates must be ascendingly sorted
    mean_south_pole_moho_depth = np.mean(depth_points[0,:])
    mean_80N_moho_depth = np.mean(depth_points[-9,:])
    mean_80S_moho_depth = np.mean(depth_points[9,:])

    grad_in_depth_N=(mean_north_pole_moho_depth - mean_80N_moho_depth)/10.0 #not sure what way round signs should go but it's consistent!
    grad_in_depth_S=(mean_south_pole_moho_depth - mean_80S_moho_depth)/-10.0

    if smooth_poles == True:

        for i in range(0,360):
            depth_points[0,i] = mean_north_pole_moho_depth
            depth_points[-1,i] = mean_south_pole_moho_depth

            for j in range(170,179): #area south of north pole
                depth_points[j,i] = mean_north_pole_moho_depth - grad_in_depth_N*(180-j)

            for k in range(0,9): #area north of south pole
                depth_points[k,i] = mean_south_pole_moho_depth - grad_in_depth_S*(k)

    if smooth_overall ==True:
    #THIS FUNCTION SMOOTHS THE WHOLE GRID: BEAR IN MIND THAT IF YOU USE THIS AFTER SMOOTHING THE POLES, IT WILL MAKE ALL THE POINTS AT
    #90 N/S HAVE SLIGHTLY DIFFERENT VALUES DUE TO THE TAILS OF THE GAUSSIAN

        depth_points = gaussian_filter(depth_points,sigma)


    if resample==True:
        #JUST IN CASE RESAMPLING THE MODEL TO HIGHER RESOLUTION REMOVED THE INSTABILITY (DIDN'T SEEM TO?)
        depth_points = zoom(depth_points, resample_factor)
        latitude_points = np.linspace(-90,90,num=180*resample_factor)
        longitude_points = np.linspace(-180, 180,num=360*resample_factor)


    try:
        ncfile = Dataset(write_file,'w',format='NETCDF4_CLASSIC')
    except PermissionError:
        Dataset(write_file).close()
        ncfile = Dataset(write_file,'w',format='NETCDF4_CLASSIC')

    ncfile.createDimension('latitude',180*resample_factor) #creates a dimension within the netcdf file with length 180*spd
    ncfile.createDimension('longitude',360*resample_factor)
    ncfile.createDimension('relative_moho_radius',0) #means that this dimension can grow as needed.

    latitude = ncfile.createVariable('latitude','d',('latitude',)) #creates variables of type 'float' within the netcdf file, with lengths and headers
    longitude = ncfile.createVariable('longitude','d',('longitude',))
    relative_moho_radius = ncfile.createVariable('relative_moho_radius','d',('latitude','longitude'))

    longitude[:] = longitude_points #actually fill data in
    latitude[:] = latitude_points

    radius_points = 6371000 + depth_points #depth points are negative

    relative_moho_radius[:] = radius_points - 6336000 ##### moho depth in the velocity model!!!!!! was 6346600 for prem

    ncfile.close()

def read_data_from_netcdf_file(filename='moho.nc', plot=False):

    model = Dataset(filename)

    latitude = model.variables['latitude']
    longitude = model.variables['longitude']
    rmr = model.variables['relative_moho_radius']

    if plot ==True:
        plot_moho(latitude, longitude, rmr)

    return model, latitude, longitude, rmr

def plot_moho(latitude,longitude, rmr):
    plt.figure(figsize=(10, 6))
    p,t = np.meshgrid(longitude, latitude)
    # ax = plt.subplot(111, projection=ccrs.PlateCarree())
    ax = plt.subplot(111, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65)))
    mesh = ax.pcolormesh(p, t, rmr, cmap='cmc.vik_r', norm=mynorm, transform=ccrs.PlateCarree())
    # mesh = ax.pcolormesh(p, t, rmr, cmap='cmc.lapaz_r', transform=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', color='white', linewidth=.5)
    ax.add_feature(cfeature.STATES, edgecolor='gray', linewidth=0.5)
    # ax.set_extent([-166.5, -136.5, 53.5, 73.5], crs=ccrs.PlateCarree())
    # ax.set_extent([-176.5, -126.5, 48.5, 78.5], crs=ccrs.PlateCarree())
    ax.set_extent([-179.5, -131.5, 37.5, 87.5], crs=ccrs.PlateCarree()) # values taken from merge_dip_moho.py

    # ax.set_extent([-129, -70, 25, 55], crs=ccrs.PlateCarree()) # for US

    # plt.figure()
    # plt.pcolormesh(p,t,rmr,cmap='cmc.lapaz_r', norm=mynorm)
    # plt.xlim([-168, -135])  # Longitude range
    # plt.ylim([50, 75])
    # cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', label="Moho (km)",fraction=.025,pad=.05)
    cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', label="Moho perturbations from 35km (m)",fraction=.025,pad=.05)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 12, 'color': 'gray'}
    gl.ylabel_style = {'size': 12, 'color': 'gray'}
    # plt.colorbar(label="moho (km)")

def plot_moho_polar_proj(latitude,longitude, rmr):
    plt.figure(figsize=(10, 6))
    p,t = np.meshgrid(longitude, latitude)
    # ax = plt.subplot(111, projection=ccrs.PlateCarree())

    ax = plt.subplot(111, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65)))
    mesh = ax.pcolormesh(p, t, rmr, cmap='cmc.vik_r', norm=mynorm, transform=ccrs.PlateCarree())
    # mesh = ax.pcolormesh(p, t, rmr, cmap='cmc.lapaz_r', transform=ccrs.PlateCarree())

    ax.coastlines(resolution='10m', color='white', linewidth=.5)
    ax.add_feature(cfeature.STATES, edgecolor='gray', linewidth=0.5)
    # ax.set_extent([-166.5, -136.5, 53.5, 73.5], crs=ccrs.PlateCarree())
    # ax.set_extent([-176.5, -126.5, 48.5, 78.5], crs=ccrs.PlateCarree())
    ax.set_extent([-179.5, -131.5, 37.5, 87.5], crs=ccrs.PlateCarree()) # values taken from merge_dip_moho.py

    # ax.set_extent([-129, -70, 25, 55], crs=ccrs.PlateCarree()) # for US

    # plt.figure()
    # plt.pcolormesh(p,t,rmr,cmap='cmc.lapaz_r', norm=mynorm)
    # plt.xlim([-168, -135])  # Longitude range
    # plt.ylim([50, 75])
    # cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', label="Moho (km)",fraction=.025,pad=.05)
    cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', label="Moho perturbations from 35km (m)",fraction=.025,pad=.05)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 12, 'color': 'gray'}
    gl.ylabel_style = {'size': 12, 'color': 'gray'}
###


#########
# write_data_to_netcdf_file()
# write_data_to_netcdf_file(smooth_overall=False, smooth_poles=False,resample=True, sigma=1, resample_factor =10, write_file = 'moho_rf_merge_no_smooth.nc')
# write_data_to_netcdf_file(smooth_overall=True, smooth_poles=False,resample=False, sigma=1, resample_factor =10, write_file = 'MM_c1_1deg_rf_smooth.nc')

write_data_to_netcdf_file(smooth_overall=False, smooth_poles=False,resample=True, sigma=1, resample_factor = 10, write_file = 'dipping_moho/moho_dip_prll_N_no_smooth.nc')

# sys.exit()

plt.ion()
# a,b,c,d=read_data_from_netcdf_file(filename='moho_rf_merge_no_smooth.nc', plot=True)
a,b,c,d=read_data_from_netcdf_file(filename='dipping_moho/moho_dip_prll_N_no_smooth.nc', plot=True)

plt.title('Dip Moho prll N not smooth')
plt.show()
sys.exit()

# plt.savefig('c1.0_moho_tr_2deg_mesh_nc_not_smooth.png', dpi=300,bbox_inches='tight', pad_inches=0.1)
# plt.savefig('merge_moho_tr_2deg_mesh_nc_not_smooth.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

plt.savefig('dipping_moho/moho_dip_prll_N_no_smooth.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

sys.exit()

a,b,c=load_moho_data(plot=True)
plt.title('Crust 1.0')
# plt.savefig('merge_moho_tr_2deg_mesh.png', dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.savefig('c1.0_US_mesh.png', dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()
