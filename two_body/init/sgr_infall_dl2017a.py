'''
Author:	Thorsten Tepper Garcia

This is setting is an attempt to reproduce the Sgr infall model of Dierickx & Loeb (2017a),
with the following differences:

- Time integration is carried on using a leapfrog (rather than RK4) method;
- MW does not include a bulge or a stellar disc; their mass is added to the DM halo;
- The scale length of the MW halo has been reduced from 38.35 to 25.
- The mass of Sgr has been increased from 1.3e10 to 2e10 to roughly match its orbital
history

Therefore, some differences between their results and ours are expected.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 7.85														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1 = 1.25e12 #1.325e12												# total mass (Msun)
a1 = 25. #38.35														# scale radius (kpc)
Potential1 = funcs.Hernquist_Potential(mass=Mass1,a=a1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
Mass2 = 1.2e10													# total mass (Msun)
a2 = 9.81														# scale radius (kpc)
Potential2 = funcs.Hernquist_Potential(mass=Mass2,a=a2)			# potential (km/s)^2
x2_0 = 125.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = -11.6													# velocities (km/s):
vy2_0 = 0.
vz2_0 = 71.66


# Dynamical friction settings
soft_length2 = 0.3 #1.0													# softening length of Sgr
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function
Mass1_cum = funcs.Hernquist_Mass(mass=Mass1,a=a1)
Mass2_cum = funcs.Hernquist_Mass(mass=Mass2,a=a2)

# Mass loss
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function

