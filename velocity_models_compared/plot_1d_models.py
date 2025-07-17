import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

taup_iasp91 = np.genfromtxt('taup_iasp91.txt')
matplotlib.rcParams.update({'font.size': 11})
######
earth_radius = 6371.0

# Layer boundaries (radii in km)
radii = np.array([6371.0, 6351.0, 6336.0, 6251.0, 6161.0, 5961.0,
                  5711.0, 5611.0, 3631.0, 3482.0, 1217.0, 0.0])

# Convert radii to depths
# depths = earth_radius - radii
depths = radii/earth_radius


# Polynomial coefficients for Vp (by layer)
vp_coeffs = [
    [5.8,0,0,0],
    [6.5,0,0,0],
    [8.78541, -0.74953,0,0],
    [25.41389, -17.69722,0,0],
    [30.78765, -23.25415,0,0],
    [29.38896, -21.40656,0,0],
    [25.96984, -16.93412,0,0],
    [25.1486, -41.1538, 51.9932, -26.6083],
    [14.49470, -1.47089,0,0],
    [10.03904, 3.75665, -13.67046],
    [11.24094, 0.0, -4.09689,0]
]
vs_coeffs = [
    [3.36,0,0,0],
    [3.75,0,0,0],
    [6.706231, -2.248585,0,0],
    [5.75020, -1.2742,0,0],
    [15.24213, -11.08552,0,0],
    [17.70732, -13.50652,0,0],
    [20.76890, -16.53147,0,0],
    [12.9303, -21.2590, 27.8988, -14.1080],
    [8.16616, -1.58206,0,0],
    [0.0,0,0,0],
    [3.56454, 0.0, -3.45241,0],
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
depth_plot=6371-depth_plot*6371
# sys.exit()
##
fig = plt.figure(figsize=(6,8))
ax = plt.axes()

# Setting the background color of the plot
# using set_facecolor() method
# ax.set_facecolor("oldlace")
ax.set_facecolor("floralwhite")

# plt.style.use('ggplot')
#plt.plot(crst1_Vp[:,1],crst1_Vp[:,0], 'darkkhaki',label='Crust 1.0')


plt.plot(taup_iasp91[:,1],taup_iasp91[:,0], 'indianred', label='Taup')

plt.plot(taup_iasp91[:,2],taup_iasp91[:,0], 'indianred')
#plt.plot(Vs_amb_YC[:,1],Vs_amb_YC[:,0], 'steelblue',linestyle='-.',label='ANT LE')

plt.plot(vp_plot, depth_plot,'yellowgreen', label="Axisem poly 3",linestyle='--')
plt.plot(vs_plot, depth_plot,'yellowgreen',linestyle='--')

######
kw={'markersize':7,'alpha':0.85,'markeredgecolor':'black','markeredgewidth':0.25}

#plt.plot(3.78,25,'o',markerfacecolor='darkkhaki',label='6.54, 3.78 km/s',**kw)
#plt.plot(6.56,25,'o',markerfacecolor='darkkhaki',**kw)

plt.xlabel('velocity (km/s)')
plt.ylabel('depth (km)')
plt.grid(color='black', linestyle='-', linewidth=.15,alpha=.35)
plt.gca().invert_yaxis()
plt.xticks(np.arange(0, 15, 1))
# plt.yticks(np.arange(0, 50, 5))
# plt.ylim((3000,0))
# plt.xlim((2,8.5))
plt.legend(fontsize="11")
plt.savefig('iasp91_comp.png',bbox_inches='tight', pad_inches=0.1)
plt.show()
