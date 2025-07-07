import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
def circle_name(index):
    """Generate names: A, B, ..., Z, AA, AB, etc."""
    name = ""
    while True:
        index, r = divmod(index, 26)
        name = chr(65 + r) + name
        if index == 0:
            break
        index -= 1
    return name


ring_radii = list(range(4, 51, 2))  # 4° to 40° every 2°
azimuths = list(range(0, 360, 10))  # 0°, 30°, ..., 340°
network = "II"


with open("stations_new.txt", "w") as f:
    f.write("# ------------------------------------------------------------------------------\n")
    f.write("#     name   network            latitude           longitude  useless      depth\n")
    f.write("# ------------------------------------------------------------------------------\n")

    for i, radius in enumerate(ring_radii):
        ring_label = circle_name(i)
        for azi in azimuths:
            lat = radius * math.cos(math.radians(azi))
            lon = radius * math.sin(math.radians(azi))
            name = f"{ring_label}{azi}"
            f.write(f"{name:>9}  {network:<6}  {lat:>18.6f}  {lon:>18.6f}  {0.0:>8.1f}  {0.0:>8.1f}\n")


###

stations = []
with open("stations_new.txt", "r") as f:
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
ax.scatter(lons, lats, color='firebrick', s=.75,alpha=.85, transform=ccrs.PlateCarree(), label='Stations')
plt.title('Concentric Ring Stations around (0°, 0°)')
plt.legend()
plt.savefig('st_cocentric.png',dpi=300,bbox_inches='tight', pad_inches=0.1)
plt.show()
