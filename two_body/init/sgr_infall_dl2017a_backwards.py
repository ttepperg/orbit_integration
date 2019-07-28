'''
Author:	Thorsten Tepper Garcia

Experimental!

Present-day mass = observed mass or final bound mass given by integrating forward

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -7.85														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1_scale = 1.25e12 #1.325e12									# mass scale (Msun)
a1 = 25. #38.35													# scale radius (kpc)
Potential1 = funcs.Hernquist_Potential(Mass1_scale,a1)	# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
Mass2_scale = 1.2e10											# mass scale (Msun)
a2 = 9.81														# scale radius (kpc)
Potential2 = funcs.Hernquist_Potential(Mass2_scale,a2)	# potential (km/s)^2
x2_0 = 12.58921													# positions (kpc)
y2_0 = 0.
z2_0 = 24.20118
vx2_0 = -365.14238												# velocities (km/s):
vy2_0 = 0.
vz2_0 = -34.97262


# Dynamical friction settings
soft_length2 = 0.3 #1.0												# softening length of Sgr
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

# Mass gain
rtrunc2 = 1.e10														# ~ inf truncation radius
Mass2_min = 6.176917e8												# present-day mass
Mass1_cum = funcs.Hernquist_Mass(Mass1_scale,a1)
Mass2_cum = funcs.Hernquist_Mass(Mass2_scale,a2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
