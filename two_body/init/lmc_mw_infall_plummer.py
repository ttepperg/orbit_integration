'''
Author:	Thorsten Tepper Garcia
Date:	11/07/2019

Experimental!

This set of parameters is used to mimick the infall of the Large Magellanic Cloud (LMC)
onto the Milky Way (MW).

The MW is approximated by a static Plummer DM halo; the LMC by a massive point-like object.

Dynamical friction by the MW is not included.


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.02e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1 = 8.e11
rs = 1.6e1														# NFW scale radius (kpc)
Potential1 = funcs.Plummer_Potential(mass=Mass1,a=rs)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
Mass2 = 1.0e10													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = -1.														# positions (kpc)
y2_0 = -41.
z2_0 = -28.
vx2_0 = 79.														# velocities (km/s):
vy2_0 = 227.
vz2_0 = -208.

