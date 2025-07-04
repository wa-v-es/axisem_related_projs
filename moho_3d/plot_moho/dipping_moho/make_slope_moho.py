# this scripts generates a box at an angle (st it is perpedicular to a chosen ray path)
# then it adds moho values along the rotated box edge, from one side to another.

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature as cfeature
from shapely.geometry import Polygon, Point
import sys
import cmcrameri as cm
from shapely.geometry import LineString

#########

# the corner points of box
corners = np.array([
    [80, -172],
    [45, -172],
    [45, -138],
    [80, -138]
])

moho_min=-5
moho_max=-65
rotation=-45
# Function to rotate points around a center
def rotate(points, angle, center):
    angle = np.deg2rad(angle)
    rot_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    return np.dot(points - center, rot_matrix) + center

###

# Rotate the box by 45 degrees
center = np.mean(corners, axis=0)
rotated_corners = rotate(corners, rotation, center)

# Create the mesh using the rotated points of the box
lat_min, lat_max = min(rotated_corners[:, 0]), max(rotated_corners[:, 0])
lon_min, lon_max = min(rotated_corners[:, 1]), max(rotated_corners[:, 1])
latitudes = np.arange(np.floor(lat_min)+.5, np.ceil(lat_max) + .5, 1)
longitudes = np.arange(np.floor(lon_min)+.5, np.ceil(lon_max) + .5, 1)
mesh_grid = np.array(np.meshgrid(latitudes, longitudes)).T.reshape(-1, 2)

# Filter mesh points to keep only those inside the rotated box
polygon = Polygon(rotated_corners)
filtered_mesh_grid = [point for point in mesh_grid if polygon.contains(Point(point))]
filtered_mesh_grid = np.array(filtered_mesh_grid)
## add moho values to filtered mesh..
lat_min, lat_max = min(filtered_mesh_grid[:, 1]), max(filtered_mesh_grid[:, 1])

# Create a line along the length of the box
# line = LineString([rotated_corners[3], rotated_corners[0]])
# line = LineString([rotated_corners[0], rotated_corners[3]]) # change the slope direction

line = LineString([rotated_corners[1], rotated_corners[0]]) # to change the slope azimuth by 90

# Calculate the linear values for the third column along the length of the box
linear_values = np.linspace(moho_min, moho_max, num=len(filtered_mesh_grid))

# Calculate the distance of each point from the start of the line
distances = np.array([line.project(Point(point)) for point in filtered_mesh_grid])

# Normalize distances to the range of linear values
normalized_distances = (distances - distances.min()) / (distances.max() - distances.min())
third_column = normalized_distances * (moho_max + -1*moho_min) - -1*moho_min

filtered_mesh_grid_with_moho = np.column_stack((filtered_mesh_grid, third_column))


fig, ax = plt.subplots(figsize=(10, 6),subplot_kw={'projection': ccrs.PlateCarree()})
# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50)})

# fig=plt.figure(figsize=(10, 6))
# ax = plt.subplot(111, projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65)))

#projection=ccrs.AlbersEqualArea(central_longitude=-154, central_latitude=50,standard_parallels=(55,65))
ax.set_extent([-180, -120, 25, 90], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.STATES, edgecolor='gray', linewidth=0.5)
ax.coastlines(resolution='50m', linewidth=.5)

# Plot the rotated box and mesh
rotated_corners = np.vstack([rotated_corners, rotated_corners[0]])
ax.plot(rotated_corners[:, 1], rotated_corners[:, 0], 'k--',alpha=.5)

mynorm = plt.Normalize(vmin=moho_max, vmax=moho_min)

sc = ax.scatter(filtered_mesh_grid_with_moho[:,1], filtered_mesh_grid_with_moho[:,0], c=filtered_mesh_grid_with_moho[:,2], cmap='cmc.vik_r', norm=mynorm, s=25, transform=ccrs.PlateCarree())
cbar = plt.colorbar(sc,ax=ax,orientation='horizontal',fraction=.035,pad=.1)

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 12, 'color': 'gray'}
gl.ylabel_style = {'size': 12, 'color': 'gray'}

# plt.title('Moho perp W deep')
plt.title('Moho prll N deep')

plt.savefig('dip_moho_prll_N_deep.png', dpi=300,bbox_inches='tight', pad_inches=0.1)

# plt.show()
##
np.savetxt('dip_moho_prll_N_deep.txt', filtered_mesh_grid_with_moho[:, [1, 0, 2]], fmt='%.1f', header='Longitude, Latitude, Moho depth')
###
