'''
Author:	Thorsten Tepper Garcia

Experimental!

Present-day mass = observed mass or final bound mass given by integrating forward

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.135e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
a1 = 1.e0														# PITS scale radius
rho01 = 2.e8
Mass1_scale = 4. * pi * rho01 * a1**3
Potential1 = funcs.PITS_Potential(Mass1_scale,a1)						# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
a2 = 1.e0														# PITS scale radius
rho02 = 1.2e8
Mass2_scale = 4. * pi * rho02 * a2**3
Potential2 = funcs.PITS_Potential(Mass2_scale,a2)				# potential (km/s)^2
x2_0 = 0.														# positions (kpc)
y2_0 = -54.65502
z2_0 = -45.81015
vx2_0 = 0.														# velocities (km/s):
vy2_0 = -53.52208
vz2_0 = 71.72406


# Dynamical friction settings
soft_length2 = 2.0e1												# softening length of Sgr
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

# Mass loss
rtrunc1 = 1e3
rtrunc2 = 70.
Mass2_min = 4.779118e10
Mass1_cum = funcs.PITS_Mass(Mass1_scale,a1)							# mass function
Mass2_cum = funcs.PITS_Mass(Mass2_scale,a2)							# mass function
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
