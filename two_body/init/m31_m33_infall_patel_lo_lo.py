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

M31 -> low mass case (model a in that paper)
M33 -> low mass case


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -6.135e0													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
c1 = 9.56														# NFW concentration
rs1 = 31.278													# NFW scale radius (kpc)
rho01 = 2.87e6													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)			# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)									# total ('virial') mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
Mass2 = 5.e10													# total mass (Msun)
rs2 = 1.															# scale radius (kpc)
Potential2 = funcs.Plummer_Potential(mass=Mass2,a=rs2)		# potential (km/s)^2
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = -24.													# velocities (km/s):
vy2_0 = 177.
vz2_0 = 94.


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)	# dynamical friction function


# Info
print("Virial radius (M31; kpc): {:.3f}".format(Rvir1))
print("Virial Mass (M31; kpc): {:E}".format(Mass1))		# should be ~2x10^12 Msun

