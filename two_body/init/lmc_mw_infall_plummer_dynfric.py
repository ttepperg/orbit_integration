'''
Author:	Thorsten Tepper Garcia

Experimental!

This is essentially identical to lmc_mw_infall_plummer.py but includes
the effect of dynamical friction of the MW onto the LMC.


'''

from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.02e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1_scale = 1.5e12
rs = 1.6e1														# NFW scale radius (kpc)
Potential1 = funcs.Plummer_Potential(Mass1_scale,rs)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Cloud
Mass2_scale = 5.0e10											# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)				# potential (km/s)^2
x2_0 = -1.														# positions (kpc)
y2_0 = -41.
z2_0 = -28.
vx2_0 = 79.														# velocities (km/s):
vy2_0 = 227.
vz2_0 = -208.


# Dynamical friction settings
soft_length2 = 5.0													# softening length of satellite
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

