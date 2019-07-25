'''
Author:	Thorsten Tepper Garcia


This set of parameters is essentially identical to m31_m33_infall_semczuk_vandermarel.py,
but taking dynamical friction into account. See m31_m33_infall_semczuk_vandermarel.py for
more information.

The resulting orbital parameters at -10 Gyr are:


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.022e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 11.65														# NFW scale radius (kpc)
rho01 = 4.2e7													# core density (Msun/kpc**3)
Mass1 = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1,rs1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
rs2 = 17.88														# NFW scale radius (kpc)
# rho02 = 3.88e6													# core density (Msun/kpc**3)
rho02 = 4.45e6													# core density (Msun/kpc**3)
Mass2 = 4. * pi * rho02 * rs2**3
Potential2 = funcs.NFW_Potential(Mass2,rs2)			# potential (km/s)^2
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = -23.2													# velocities (km/s):
vy2_0 = 177.4
vz2_0 = 93.7


# Dynamical friction settings
soft_length2 = 28.5													# softening length of M33
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

