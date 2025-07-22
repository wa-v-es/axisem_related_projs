# originally from https://github.com/AxiSEMunity/AxiSEM3D/blob/master/examples/04_simple_3d_shapes/Example_4_README.ipynb
## The plume tail axis at [0, 0] in coordinates and a depth from 1000 km to 3000 km.
# The plume head is a sphere (ellipsoid) of radius 500 km centred at a depth of 500 km.
#The Vp, Vs and density perturbation will all be - 10 %.
import sys

import numpy as np
import scipy.ndimage as sci
sys.path.append('src/')
from model import Model
from ellipsoid import Ellipsoid
from cylinder import *
from injector import *
import netCDF4 as nc


# We dont need to make a model that spans the whole domain, just the part we are interested in injecting a plume in:
radius = 2800000
perturb = -0.1
lat_lim = [-50, 50]
long_lim = [-50, 50]
depth_lim = [0, radius]

# Set locations for shapes:
ell_loc = [0,0, 400000] # orignally was 1200km.
cyl_loc = [0,0, 2200000] # [x,y,z] of centre of cylinder
# sys.exit()
# Create our global model:
glob_m = Model("spherical", lat_lim, long_lim, depth_lim, elements_per_wavelength=1, dominant_freq=.5, min_velocity=10000, oversaturation=1, a=radius)

# Create cylinder: [h, rad, theta, phi, expand_int] where h - length of the cylinder, radius of the cylinder, theta and phi are rotation angles away from the major axis.
cylinder = Cylinder(model=glob_m, vp=perturb, vs=perturb, rho=perturb, dim=[2500000, 300000, 0, 0, 1], loc=cyl_loc, major_axis='Z')

# Create ellipse:
ellipse = Ellipsoid(model=glob_m, vp=perturb, vs=perturb, rho=perturb, dim=[500000, 500000, 500000, np.pi/2, 0, 1], loc=ell_loc)

# Create injector object and inject
i = Injector(glob_m)
i.addObj(cylinder, location=cyl_loc, overwrite=True)
i.addObj(ellipse, location=ell_loc, overwrite=True)

# Gaussian filter to make the plume boundaries a bit less harsh
sigma = 1
glob_m.bm_rho =  sci.gaussian_filter(input=glob_m.bm_rho, sigma=sigma)
glob_m.bm_vp =  sci.gaussian_filter(input=glob_m.bm_vp, sigma=sigma)
glob_m.bm_vs =  sci.gaussian_filter(input=glob_m.bm_vs, sigma=sigma)

# Write to NetCDF file
out_path = f"plume_{sigma}"
glob_m.writeNetCDF(f"{out_path}mantle.nc")


glob_m.writeNetCDF(f"{out_path}mantle_visual.nc", paraview=True)

###
