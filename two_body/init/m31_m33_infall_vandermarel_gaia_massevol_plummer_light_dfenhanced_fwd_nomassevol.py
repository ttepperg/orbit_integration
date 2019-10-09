'''
Author:	Thorsten Tepper Garcia

This setup is essentially identica to
m31_m33_infall_vandermarel_gaia_massevol_plummer_light_dfenhanced_fwd.py
but ignoring tidal stripping.


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
rs2 = 17.88														# Plummer scale radius (kpc)
rho02 = 3.71e6													# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)	# potential (km/s)^2
x2_0 = -320.71800													# positions (kpc)
y2_0 = 306.05750
z2_0 = -172.37540
vx2_0 = 174.80037													# velocities (km/s):
vy2_0 = -107.14605
vz2_0 = 115.53696

# Dynamical friction settings
# eps is tuned to match Nbody sim
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=2.e-1*rs2)
