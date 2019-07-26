'''
Author:	Thorsten Tepper Garcia

This setup is an attempt to model the infall of the Sagittarius
dwarf galaxy onto the Milky Way. The observed data used as initial conditons
correspond to the values quoted by Dierickx & Loeb (2017a).
We also adopt their mass models for each galaxy.

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
x2_0 = 16.1														# positions (kpc)
y2_0 = 2.35
z2_0 = -6.12
vx2_0 = 242.5													# velocities (km/s):
vy2_0 = 5.6
vz2_0 = 228.1


# Dynamical friction settings
soft_length2 = 1.0													# softening length of LMC
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

