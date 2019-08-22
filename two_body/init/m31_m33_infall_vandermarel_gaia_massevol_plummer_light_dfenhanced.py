'''
Author:	Thorsten Tepper Garcia

This setup is essentially identical to
m31_m33_infall_vandermarel_gaia_massevol_plummer_light_fwd.py, but
adopting a different (smaller) smoothing length for M33.
The reason is that an N-body simulation based on the former suggests
that the semi-analytic calculation severely underestimates the
effect of dynamical friction OR overestimates tidal stripping.


The infall state vectors are:

	Infall relative position: (-320.71800,306.05750,-172.37540)
	Infall relative distance: 475.65167
	Infall relative velocity: (174.80037,-107.14605,115.53696)
	Infall relative speed: 235.33855

	Infall bound mass of body 1: 2.003856E+12
	Infall bound mass of body 2: 2.306989E+11

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -8.e0														# total time (time unit ~ 0.978 Gyr)
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
rs2 = 17.88														# Plummer scale radius (kpc)
rho02 = 5e6														# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)				# potential (km/s)^2
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = 49.														# velocities (km/s):
vy2_0 = 190.
vz2_0 = 112.


# Dynamical friction settings
# eps is tuned to match Nbody sim
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=2.e-1*rs2)

# Mass evolution (unstripping)
Mass2_min = 1.3e11															# present-day mass
Mass1_cum = funcs.NFW_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Plummer_Mass(Mass2_scale,rs2)
Mass2_evol = funcs.mass_bound(Mass1_cum,Mass2_cum)					# mass evolution function
