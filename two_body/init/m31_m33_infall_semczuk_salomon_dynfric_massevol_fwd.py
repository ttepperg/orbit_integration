'''
Author:	Thorsten Tepper Garcia

This setup is essentially aims at tuning the infall of M33 onto M31
iteratively by hand, alternating between this setup and
m31_m33_infall_semczuk_salomon_dynfric_massevol.py.


The boundary conditions are:
- present-day position and velocity of M33 relative to M31
- M33's present-day mass => ~3x10^11 Msun (at the low end of the range
  estimated by Corbelli et al. 2014)


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 5.112														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 11.65														# NFW scale radius (kpc)
rho01 = 4.2e7													# core density (Msun/kpc**3)
rtrunc1 = 326.													# truncation (virial) radius
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)				# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
rs2 = 17.88														# NFW scale radius (kpc)
rho02 = 1e7														# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.NFW_Potential(Mass2_scale,rs2)				# potential (km/s)^2
x2_0 = 56.19431													# positions (kpc)
y2_0 = -70.18257
z2_0 = -9.93614
vx2_0 = -329.37955												# velocities (km/s):
vy2_0 = -107.56652
vz2_0 = -255.66212


# Dynamical friction settings
soft_length2 = 28.5													# softening length of M33
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

# Mass evolution (stripping)
Mass1_cum = funcs.Hernquist_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Hernquist_Mass(Mass2_scale,rs2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
