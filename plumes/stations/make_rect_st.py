import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#
network = "II"
lat_range = range(-20, 21, 1)  # -20° to 20°
lon_range = range(5, 21, 1)    # 4° to 30°

def format_station_name(lat, lon):
    lat_dir = 'N' if lat >= 0 else 'S'
    lon_dir = 'E'  # all positive in your range
    return f"la{abs(lat)}{lat_dir}lo{abs(lon)}{lon_dir}"

with open("grid_stations.txt", "w") as f:
    f.write("# ------------------------------------------------------------------------------\n")
    f.write("#     name   network            latitude           longitude  useless      depth\n")
    f.write("# ------------------------------------------------------------------------------\n")

    for lat in lat_range:
        for lon in lon_range:
            name = format_station_name(lat, lon)
            f.write(f"{name:>9}  {network:<6}  {lat:>18.6f}  {lon:>18.6f}  {0.0:>8.1f}  {0.0:>8.1f}\n")

print("Station data has been written to 'grid_stations.txt'.")

stations = []

with open("grid_stations.txt", "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        lat = float(parts[2])
        lon = float(parts[3])
        stations.append((lat, lon))

# Extract latitudes and longitudes
lats, lons = zip(*stations)

# Plotting
plt.figure(figsize=(3, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.stock_img()
ax.coastlines()
ax.gridlines(draw_labels=True)
ax.scatter(lons, lats, color='firebrick', s=.75,alpha=.85, transform=ccrs.PlateCarree(), label='Stations')
plt.title('Grid Stations')
plt.legend()
plt.savefig('st_grid.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()
