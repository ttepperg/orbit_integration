'''
Author:	Thorsten Tepper Garcia

This set of parameters is essentially identical to m31_m33_infall_semczuk_salomon.py,
but taking dynamical friction into account. See m31_m33_infall_semczuk_salomon.py for
more information.

The resulting orbital parameters at -5 Gyr are:


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -5.112													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
c1 = 28.														# NFW concentration
rs1 = 11.65														# NFW scale radius (kpc)
rho01 = 4.2e7													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)						# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)										# total ('virial') mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
c2 = 11.														# NFW concentration
rs2 = 17.88														# NFW scale radius (kpc)
# rho02 = 3.88e6													# core density (Msun/kpc**3)
rho02 = 4.45e6													# core density (Msun/kpc**3)
Rvir2 = c2*rs2
Potential2 = funcs.NFW_Potential(rho02,rs2)						# potential (km/s)^2
Mass2_cum = funcs.NFW_Mass(rho02,rs2)
Mass2 = Mass2_cum(Rvir2)										# total ('virial') mass (Msun)
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = -72.0													# velocities (km/s):
vy2_0 = 86.4
vz2_0 = 10.6


# Dynamical friction settings
soft_length2 = 28.5													# softening length of M33
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Info
print("Virial radius (M31; kpc): {}".format(Rvir1))
print("Virial Mass (M31; kpc): {:E}".format(Mass1))		# should be ~2x10^12 Msun
print("Virial radius (M33; kpc): {}".format(Rvir2))
print("Virial Mass (M33; kpc): {:E}".format(Mass2))		# should be ~5x10^11 Msun

# exit()