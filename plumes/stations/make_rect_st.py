import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.patches import Circle
#
network = "II"
lat_range = range(-15, 16, 1)  # -20° to 20°
lon_range = range(5, 21, 1)    # 4° to 30°

def format_station_name(lat, lon):
    lat_dir = 'N' if lat >= 0 else 'S'
    lon_dir = 'E'  # all positive in your range
    return f"la{abs(lat)}{lat_dir}lo{abs(lon)}{lon_dir}"

with open("grid_stations_small.txt", "w") as f:
    f.write("# ------------------------------------------------------------------------------\n")
    f.write("#     name   network            latitude           longitude  useless      depth\n")
    f.write("# ------------------------------------------------------------------------------\n")

    for lat in lat_range:
        for lon in lon_range:
            name = format_station_name(lat, lon)
            f.write(f"{name:>9}  {network:<6}  {lat:>18.6f}  {lon:>18.6f}  {0.0:>8.1f}  {0.0:>8.1f}\n")

print("Station data has been written to 'grid_stations.txt'.")

stations = []

with open("grid_stations_small.txt", "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        lat = float(parts[2])
        lon = float(parts[3])
        stations.append((lat, lon))

# Extract latitudes and longitudes
lats, lons = zip(*stations)

circle = Circle((0, 0),
    radius=4.5,          # degree
    transform=ccrs.PlateCarree(),facecolor='khaki',edgecolor='grey',linewidth=.1,alpha=0.85,
    label='Plume'
)
cylinder = Circle((0, 0),
    radius=2.6,          #deg
    transform=ccrs.PlateCarree(),facecolor='darkkhaki',edgecolor='grey',linewidth=.1,alpha=0.5
)
# Plotting
plt.rcParams['font.size'] = 13

plt.figure(figsize=(7, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.stock_img()
ax.coastlines()
# ax.gridlines(draw_labels=True,zorder=1)
gl = ax.gridlines(draw_labels=True,
                  linewidth=0.75, color='gray', alpha=0.65, linestyle='--')

gl.top_labels = False     # remove top longitude labels
gl.right_labels = False   # remove right latitude labels
# ax.scatter(0, 0, color='khaki',marker='o', s=2500,alpha=.85, transform=ccrs.PlateCarree(), label='plume')
ax.add_patch(circle)
ax.add_patch(cylinder)


ax.scatter(lons, lats, color='orchid',marker='^', s=15,alpha=.85, transform=ccrs.PlateCarree(), label='Stations',zorder=10)
ax.scatter(5, 5, color='forestgreen',marker='^', s=50,alpha=.85, transform=ccrs.PlateCarree(),zorder=10)
ax.scatter(15, 5, color='forestgreen',marker='^', s=50,alpha=.85, transform=ccrs.PlateCarree(),zorder=10)
ax.scatter(15, 15, color='forestgreen',marker='^', s=50,alpha=.85, transform=ccrs.PlateCarree(),zorder=10)


ax.scatter(-20, -2, color='firebrick',marker='*', s=150,alpha=.85, transform=ccrs.PlateCarree(), label='Eq',zorder=10)
#
# ax.set_xticks(range(-25, 25, 5), crs=ccrs.PlateCarree())  # bottom only
# ax.set_yticks(range(-20, 20, 5), crs=ccrs.PlateCarree())    # left only
# ax.xaxis.set_tick_params(top=False, labeltop=False)
# ax.yaxis.set_tick_params(right=False, labelright=False)

# plt.title('Grid Stations')
plt.legend()
# plt.savefig('st_plume_eq.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()
