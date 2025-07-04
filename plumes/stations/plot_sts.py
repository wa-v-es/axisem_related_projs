
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Load station data
stations = []
with open("stations_output.txt", "r") as f:
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
plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.stock_img()
ax.coastlines()
ax.gridlines(draw_labels=True)
ax.scatter(lons, lats, color='firebrick', s=.75,alpha=.75, transform=ccrs.PlateCarree(), label='Stations')
plt.title('Concentric Ring Stations around (0°, 0°)')
plt.legend()
plt.show()
