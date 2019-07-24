'''
Author:	Thorsten Tepper Garcia

This setup is essentially identical to sgr_infall_backwards.py, but taking
the mass evolution (gain) of the LMC when integrating backwards.
Recall: Because mass gets stripped at infall, the mass must have been higher
at infall. This fact is taken into account here.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1 = 1.325e12												# total mass (Msun)
a1 = 38.35														# scale radius (kpc)
Potential1 = funcs.Hernquist_Potential(mass=Mass1,a=a1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
Mass2 = 1.9e9													# total mass (Msun)
a2 = 9.81														# scale radius (kpc)
Potential2 = funcs.Hernquist_Potential(mass=Mass2,a=a2)			# potential (km/s)^2
x2_0 = 16.1														# positions (kpc)
y2_0 = 2.35
z2_0 = -6.12
vx2_0 = 242.5													# velocities (km/s):
vy2_0 = 5.6
vz2_0 = 228.1


# Dynamical friction settings
soft_length2 = 1.0													# softening length of LMC
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Mass evolution
Mass1_cum = funcs.Hernquist_Mass(mass=Mass1,a=a1)
Mass2_cum = funcs.Hernquist_Mass(mass=Mass2,a=a2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function


# Test
# t=0
# r0 = [x2_0,y2_0,z2_0]
# rt = funcs.tidal_radius(*r0,m1_func=Mass1_cum,m2_func=Mass2_cum)
# print(funcs.norm(*r0))
# print(rt)
# print("{:E}".format(Mass2_cum(rt)))
# print("{:E}".format(Mass2_evol(t,*r0)))
# exit()



