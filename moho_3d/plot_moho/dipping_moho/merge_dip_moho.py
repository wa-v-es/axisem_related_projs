# using two Moho (local and global) text files of same dimensions (lat long moho depth), creates a merged Moho file such that
# local values are taken when present, and a linear average is taken for points at the edge of local data
# (transition width controls this).

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

##
# moho from miller moresi
# moho_model = alaskamoho.MohoModel_opt
# moho_data_filename="AlaskaMohoOpt"

# map_extent=[-166.5, -136.5, 53.5, 73.5]
map_extent=[-179.5, -131.5, 37.5, 87.5]
# these values are taken from lons_lats_dip[:,0].min()/max; lons_lats_dip[:,1].min()/max
#except added changed -178.5 to -179.5.
#
# lons = np.linspace(-165.5,-135.5,100)
# lats = np.linspace(55, 72, 100)
lons=np.arange(map_extent[0],map_extent[1],1)
lats=np.arange(map_extent[2],map_extent[3],1)

# Define the map projection
data_dip = np.loadtxt("dip_moho_prll_N_deep.txt",skiprows=1)  # Replace with your actual filename
lons_lats_dip = data_dip[:, :2]
# lats = data[:, 1]
depths_dip = data_dip[:, 2]

data_c1 = np.loadtxt("depthtomoho.xyz")  # Replace with your actual filename

data_c1[:,2]=-35 # makes the global moho 35km
# sys.exit()
##
data_merged=data_c1.copy()
##########################
transition_width = 2 # 2 degree
##########################

# Identify matching global data points within local extent
global_mask = (data_c1[:, 0] >= map_extent[0]) & (data_c1[:, 0] <= map_extent[1]) & \
              (data_c1[:, 1] >= map_extent[2]) & (data_c1[:, 1] <= map_extent[3])

global_subset = data_c1[global_mask]
global_indices = np.where(global_mask)[0]

# Build a KDTree for fast nearest-neighbor searching
tree = cKDTree(data_dip[:, :2])

def blending_weight(dist, transition_width):
    return np.clip(dist / transition_width, 0, 1)  # Linear blending

# Update global data smoothly
for i, (lon, lat, g_val) in enumerate(global_subset):
    dist, idx = tree.query([lon, lat], k=1)  # Find nearest local point
    if dist < transition_width:  # Apply blending if within transition zone
        print('dist = {:.2f}'.format(dist))
        local_val = data_dip[idx, 2]
        w = blending_weight(dist, transition_width)
        # Get index in full data_merged array
        idx_global = global_indices[i]  # Map i to actual index in data_merged
        data_ori = data_merged[idx_global, 2]  # Debugging step

        # Proper in-place modification
        data_merged[idx_global, 2] = w * g_val + (1 - w) * local_val

        # data_merged[global_mask][i, 2] = w * g_val + (1 - w) * local_val

        print('w={:.2f}; data_ori={:.2f}; data_merge={:.2f}\n'.format(w,data_ori,data_merged[global_mask][i, 2]))
        print('###############\n')

# sys.exit()
##
plt.ion()
fig = plt.figure(figsize=(15, 6))
ax = plt.subplot(131, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65)))

# Add map features
ax.set_extent(map_extent, crs=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='none',edgecolor='black')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=":")

# Scatter plot depth values as circles
mynorm = plt.Normalize(vmin=-65, vmax=-5)
sm = plt.cm.ScalarMappable(cmap='cmc.lapaz_r', norm=mynorm)
# color_bar = plt.colorbar(sm)
sc = ax.scatter(lons_lats_dip[:,0], lons_lats_dip[:,1], c=depths_dip, cmap='cmc.lapaz_r', norm=mynorm, s=25, transform=ccrs.PlateCarree())
ax.set_title('Dip Moho prll N ({:.1f}, {:.1f})'.format(np.min(depths_dip),np.max(depths_dip)))

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 12, 'color': 'gray'}
gl.ylabel_style = {'size': 12, 'color': 'gray'}

##
###
### USING MERGED DATA
lons_lats_c1 = data_merged[:, :2]
depths_c1 = data_merged[:, 2]

ax1 = plt.subplot(132, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65)))
ax1.set_extent(map_extent, crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, facecolor='none',edgecolor='black')
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS, linestyle=":")
sc = ax1.scatter(lons_lats_c1[:,0], lons_lats_c1[:,1], c=depths_c1, cmap='cmc.lapaz_r', norm=mynorm, s=25, transform=ccrs.PlateCarree())

###plot difference
# Create a dictionary from the with (lat, lon) as the key and third column as the value
data_c1_dict = { (lat, lon): value for lat, lon, value in data_merged }

# Iterate through the MM moho and compute differences where lat/lon matches
differences = []
for lat, lon, value in data_dip[:,:3]:
    key = (lat, lon)
    if key in data_c1_dict:
        diff = value - data_c1_dict[key]  # Compute the difference
        differences.append((lat, lon, -1*diff,value,data_c1_dict[key]))

# Convert to NumPy array for further processing
differences = np.array(differences)
###
diff_merg_c1= data_merged[:,2] - data_c1[:,2]

# sys.exit()
ax2 = plt.subplot(133, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65)))
ax2.set_extent(map_extent, crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, facecolor='none',edgecolor='black')
ax2.add_feature(cfeature.COASTLINE)
ax2.add_feature(cfeature.BORDERS, linestyle=":")

mynorm = plt.Normalize(vmin=-30, vmax=30)
sm_diff = plt.cm.ScalarMappable(cmap='cmc.vik', norm=mynorm)

sc = ax2.scatter(lons_lats_c1[:,0], lons_lats_c1[:,1], c=-data_merged[:,2] + data_c1[:,2], cmap='cmc.vik', norm=mynorm, s=30, transform=ccrs.PlateCarree())

# ax2.set_title("MM - Crust 1.0 (km)")
ax2.set_title('Merge - C1.0 ({:.1f}, {:.1f})'.format(-1*np.min(diff_merg_c1),np.max(diff_merg_c1)))
ax1.set_title('Dip Moho + C1.0 ({:.1f}, {:.1f})'.format(np.min(differences[:,4]),np.max(differences[:,4])))

###
ax3=fig.add_axes([0.25, .1, .3, .025]) #[left, bottom, width, height]
cbar = plt.colorbar(sm,cax=ax3,orientation='horizontal')
cbar.set_label('Moho (km)', rotation=0, labelpad=2,loc='center')
###
ax4=fig.add_axes([0.71, .1, .15, .025]) #[left, bottom, width, height]
cbar = plt.colorbar(sm_diff,cax=ax4,orientation='horizontal')
cbar.set_label('Moho diff (km)', rotation=0, labelpad=2,loc='center')

np.savetxt("dip_prll_N_c1_tr2.XYZ", data_merged, fmt='%.2f %.2f %.2f')#,header="Longitude Latitude Depth quality ")

fig.savefig('dip_prll_N_c1_tr2.png', dpi=400,bbox_inches='tight', pad_inches=0.1)
plt.show()

sys.exit()
###
