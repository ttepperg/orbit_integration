'''
Author:	Thorsten Tepper Garcia

This setup is essentially m31_m33_infall_vandermarel_gaia_massevol_plummer_fwd.py,
but integrated for a shorter time span. The idea is to use the final state vectors
as input for an N-body + hydrodynamical model, and thus shorten the integration
time relative to the full orbit m31_m33_infall_vandermarel_gaia_massevol_plummer_fwd.py.

The validity of this approach is tested with m31_m33_infall_vandermarel_gaia_massevol_plummer_fwd_short1.py


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 5e0														# total time (time unit ~ 0.978 Gyr)
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
rs2 = 17.88														# Plummer scale radius (kpc)
rho02 = 6.135e6												# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)				# potential (km/s)^2
x2_0 = -944.23627												# positions (kpc)
y2_0 = 444.40255
z2_0 = -672.72794
vx2_0 = 81.41245												# velocities (km/s):
vy2_0 = -13.60903
vz2_0 = 66.94253


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)		# dynamical friction function

# Mass evolution (stripping)
Mass1_cum = funcs.NFW_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Plummer_Mass(Mass2_scale,rs2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
