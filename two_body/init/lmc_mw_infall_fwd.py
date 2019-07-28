'''
Author:	Thorsten Tepper Garcia

This setup is the forward integrated version of lmc_mw_infall.py.

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0													# initial time (Gyr)
t_1 = 1.0226e1												# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3											# integration time step

# Milky Way
rs1 = 1.6e1													# NFW scale radius (kpc)
rho01 = 1.307e7												# core density (Msun/kpc**3)
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)			# potential (km/s)^2
x1_0 = 0.													# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.													# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
a2=3.0
Mass2_scale = 5.0e10										# mass scale (Msun)
Potential2 = funcs.Plummer_Potential(Mass2_scale,a2)		# potential (km/s)^2
x2_0 = -33.69621											# positions (kpc)
y2_0 = -39.87210
z2_0 = -384.92628
vx2_0 = 16.66207											# velocities (km/s):
vy2_0 = -46.33569
vz2_0 = 139.19027


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=a2)	# dynamical friction function

