'''
Author:	Thorsten Tepper Garcia

This setup is intended to check the validity of the setup
m31_m33_infall_vandermarel_gaia_massevol_plummer_fwd_short1.py
as a short-integrated version of
m31_m33_infall_vandermarel_gaia_massevol_plummer_fwd.py


The present-day state vectors in the latter are:

	Present-day relative position: (-95.74430,-122.02125,-128.66744)
	Present-day relative distance: 201.52287
	Present-day relative velocity: (45.85365,188.81194,108.79276)
	Present-day relative speed: 222.68446

while this setup yields:

	Present-day relative position: (-95.74434,-122.02140,-128.66753)
	Present-day relative distance: 201.52304
	Present-day relative velocity: (45.85357,188.81185,108.79266)
	Present-day relative speed: 222.68432

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 8.3e0														# total time (time unit ~ 0.978 Gyr)
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
rho02 = 6.135e6												# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)				# potential (km/s)^2
x2_0 = -416.88327												# positions (kpc)
y2_0 = 307.50793
z2_0 = -256.74009
vx2_0 = 142.36660												# velocities (km/s):
vy2_0 = -51.83579
vz2_0 = 106.91843


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)		# dynamical friction function

# Mass evolution (stripping)
Mass1_cum = funcs.NFW_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Plummer_Mass(Mass2_scale,rs2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
