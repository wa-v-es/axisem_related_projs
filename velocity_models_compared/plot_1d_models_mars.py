import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

taup_iasp91 = np.genfromtxt('taup_iasp91.txt')
matplotlib.rcParams.update({'font.size': 11})
######
earth_radius = 6371.0
mars_radius = 3390.0


# Layer boundaries (radii in km)
radii = np.array([3390, 3280, 3055.5, 2908, 2360, 2033, 1468, 515, 0])
# 3390 3280 3055.5 2908 2360 2033 1468 515 0
# Convert radii to depths
# depths = earth_radius - radii
depths = radii/mars_radius


# Polynomial coefficients for Vp (by layer)

vp_coeffs = [[50.537, -123.68, 123.81, -43.012],
    [14.372, -13.395, 13.353, -6.3398],
    [15.559, -17.460, 18.004, -8.1157],
    [11.9300, -4.8826, 3.5016, -2.5538],
    [11.8365, -4.1713, 3.1366, -2.5691],
    [11.4166, -0.90421, -2.6380, 0.94287],
    [6.5395, -0.0225299, -2.3767, -0.72716],
    [6.5395, -0.0225299, -2.3767, -0.72716]
]


vs_coeffs = [
    [-43.819, 148.53, -151.21, 50.526],
    [8.6113, -11.932, 14.301, -6.5441],
    [22.395, -57.011, 63.447, -24.406],
    [6.5847, -2.5543, 1.6024, -1.2675],
    [6.5172, -1.8055, 0.8080, -0.95676],
    [6.7644, -2.3139, 1.7809, -1.5312],
    [0.0, 0, 0, 0],
    [2.4, 0, 0, 0]
]

# Prepare depth and velocity arrays for plotting
depth_plot = []
vp_plot = []
vs_plot = []


vp_min, vp_max = 0, 14

# Loop through each layer
for i in range(len(depths) - 1):
    # Define the depth range for this layer
    layer_depths = np.linspace(depths[i], depths[i + 1], 100)
    # Evaluate the polynomial for this layer
    coeffs = vp_coeffs[i]
    coeffs_vs = vs_coeffs[i]

    layer_vp = np.polyval(coeffs[::-1], layer_depths)
    layer_vs = np.polyval(coeffs_vs[::-1], layer_depths)

    # layer_vp[layer_vp < 0] = np.nan
    # layer_vp = np.clip(layer_vp, vp_min, vp_max)
    # Append to plot arrays
    depth_plot.extend(layer_depths)
    vp_plot.extend(layer_vp)
    vs_plot.extend(layer_vs)


# Convert to arrays
depth_plot = np.array(depth_plot)
vp_plot = np.array(vp_plot)
vs_plot = np.array(vs_plot)


print(f"Vp range: min={vp_plot.min()}, max={vp_plot.max()}")
## converting depth back to depth from surface..
depth_plot=3390-depth_plot*3390
# sys.exit()
##
fig = plt.figure(figsize=(6,8))
ax = plt.axes()

# Setting the background color of the plot
# using set_facecolor() method
ax.set_facecolor("azure")
# ax.set_facecolor("floralwhite")

# plt.style.use('ggplot')
#plt.plot(crst1_Vp[:,1],crst1_Vp[:,0], 'darkkhaki',label='Crust 1.0')


# plt.plot(taup_iasp91[:,1],taup_iasp91[:,0], 'indianred', label='Taup')
# plt.plot(taup_iasp91[:,2],taup_iasp91[:,0], 'indianred')
#plt.plot(Vs_amb_YC[:,1],Vs_amb_YC[:,0], 'steelblue',linestyle='-.',label='ANT LE')

plt.plot(vp_plot, depth_plot,'indianred', label="Vp",linestyle='-')
plt.plot(vs_plot, depth_plot,'yellowgreen',label="Vs",linestyle='-')

######
kw={'markersize':7,'alpha':0.85,'markeredgecolor':'black','markeredgewidth':0.25}

#plt.plot(3.78,25,'o',markerfacecolor='darkkhaki',label='6.54, 3.78 km/s',**kw)
#plt.plot(6.56,25,'o',markerfacecolor='darkkhaki',**kw)

plt.xlabel('velocity (km/s)')
plt.ylabel('depth (km)')
plt.grid(color='black', linestyle='-', linewidth=.15,alpha=.35)
plt.gca().invert_yaxis()
plt.xticks(np.arange(0, 12, 1))
# plt.yticks(np.arange(0, 50, 5))
# plt.ylim((3000,0))
# plt.xlim((2,8.5))
plt.legend(fontsize="11")
plt.savefig('mars_sohl.png',bbox_inches='tight', pad_inches=0.1)
plt.show()
