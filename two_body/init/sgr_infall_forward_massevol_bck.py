'''
Author:	Thorsten Tepper Garcia

This setup is essentially identical to sgr_infall_forward_massevol.py,
but integrating backwards in time.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1_scale = 1.325e12												# mass scale (Msun)
a1 = 38.35														# scale radius (kpc)
Potential1 = funcs.Hernquist_Potential(Mass1_scale,a1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
Mass2_scale = 1.3e10													# mass scale (Msun)
a2 = 9.81														# scale radius (kpc)
Potential2 = funcs.Hernquist_Potential(Mass2_scale,a2)			# potential (km/s)^2
x2_0 = 38.00363													# positions (kpc)
y2_0 = 20.98284
z2_0 = -180.36647
vx2_0 = 80.22624												# velocities (km/s):
vy2_0 = 20.21228
vz2_0 = -121.88746


# Dynamical friction settings
soft_length2 = 1.0													# softening length of LMC
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Mass loss
rtrunc2 = 1.e10														# ~ inf truncation radius
Mass2_min = 1.669875e09
Mass1_cum = funcs.Hernquist_Mass(Mass1_scale,a1)
Mass2_cum = funcs.Hernquist_Mass(Mass2_scale,a2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
