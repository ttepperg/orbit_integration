'''
Author:	Thorsten Tepper Garcia

This set of parameters is intended to reproduce the semi-analytic calculation
for the infall of M33 onto M31 by Patel et al. (2017) with the following differences:

- We approximate M31 by a DM halo only, and ignore additional
component (e.g. stellar bulge, stellar disc). The mass of these components
is added in each case to the host DM halo.
- The expression for the Couloumb logarithm differs slightly, more specifically,
the constant factor in the denominator we use is 1.4 rather than 1.2


The mass models for the galaxies are based on the extreme cases explored
in Patel et al. (2017). Here, we adopt:

M31 -> high mass case (model b in that paper)
M33 -> high mass case


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -6.135e0													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 35.15														# NFW scale radius (kpc)
rho01 = 2.683e6												# core density (Msun/kpc**3)
Mass1 = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1,rs1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
Mass2 = 2.5e11													# total mass (Msun)
rs2 = 21.														# scale radius (kpc)
Potential2 = funcs.Plummer_Potential(Mass2,rs2)		# potential (km/s)^2
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = -24.													# velocities (km/s):
vy2_0 = 177.
vz2_0 = 94.


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)	# dynamical friction function

