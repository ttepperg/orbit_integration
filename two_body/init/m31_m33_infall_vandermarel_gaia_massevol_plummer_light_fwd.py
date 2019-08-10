'''
Author:	Thorsten Tepper Garcia

This setup is essentially aims at tuning the infall of M33 onto M31
iteratively by hand, alternating between this setup and
m31_m33_infall_vandermarel_gaia_massevol.py.

The boundary conditions are:
- present-day position and velocity of M33 relative to M31
- M33's present-day mass => ~1.3x10^11 Msun (at the centre of the range
adopted by Patel et al. 2017).
- Infall mass of M33 should be ~2.6x10^11 Msun

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 8.e0														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 11.65														# NFW scale radius (kpc)
rho01 = 4.2e7													# core density (Msun/kpc**3)
rtrunc1 = 326.													# truncation (virial) radius
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)	# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
rs2 = 17.88														# NFW scale radius (kpc)
rho02 = 3.66e6 #7e6														# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)	# potential (km/s)^2
x2_0 = -322.98444													# positions (kpc)
y2_0 = 247.08110
z2_0 = -195.71482
vx2_0 = 121.60311													# velocities (km/s):
vy2_0 = -42.69608
vz2_0 = 91.89645


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)		# dynamical friction function

# Mass evolution (stripping)
Mass1_cum = funcs.NFW_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Plummer_Mass(Mass2_scale,rs2)
Mass2_evol = funcs.mass_bound(Mass1_cum,Mass2_cum)	# mass evolution function
