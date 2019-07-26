'''
Author:	Thorsten Tepper Garcia

Experimental!

This set of parameters is used to mimick the infall of the Large Magellanic Cloud (LMC)
onto the Milky Way (MW).

The MW is approximated by a static Hernquist DM halo; the LMC by a massive point-like object.

Dynamical friction by the MW is not included.


'''

from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.02e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1_scale = 1.5e12
rs = 4.0e1														# NFW scale radius (kpc)
Potential1 = funcs.Hernquist_Potential(Mass1_scale,rs)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
Mass2_scale = 5.0e10													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)					# potential (km/s)^2
x2_0 = -1.														# positions (kpc)
y2_0 = -41.
z2_0 = -28.
vx2_0 = 79.														# velocities (km/s):
vy2_0 = 227.
vz2_0 = -208.

