'''
Author:	Thorsten Tepper Garcia

This setting essentially corresponds to lmc_mw_infall_galpy_dynfric.py,
but integrated forward in time, using the state vectors provided by
the backwards integration as initial conditions.

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
rs1 = 1.6e1														# NFW scale radius (kpc)
rho01 = 1.307e7												# core density (Msun/kpc**3)
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
Mass2_scale = 5.0e10													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)			# potential (km/s)^2
x2_0 = -44.48375												# positions (kpc)
y2_0 = 13.77154
z2_0 = -456.73296
vx2_0 = 12.91274												# velocities (km/s):
vy2_0 = -52.52823
vz2_0 = 95.00023

# Dynamical friction settings
soft_length2 = 5.0																	# softening length of LMC
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function
