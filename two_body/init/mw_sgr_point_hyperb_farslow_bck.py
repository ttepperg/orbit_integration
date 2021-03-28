'''
Author: Thorsten Tepper Garcia
Date: 29/04/2020

Modifications
17/05/2020	- changed satellite mass from 5.6e10 to 2e10
			- changed name from mw_sgr_hyperb_bck_point_farslow.py to
				mw_sgr_point_hyperb_farslow_bck.py
			- swapped y2_0 and z2_0 coordinates
			- changed rho01 from 1.57e7 to 1.64e7


THIS is the orbit actually used to model the MW-Sgr interaction with AGAMA paper (Bland-Hawthorn & Tepper-Garcia 2021)

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0             # initial time (Gyr)
t_1 = -1.0e-1           # total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-4        # integration time step

# Milky Way
rs1 = 1.5e1             # NFW scale radius (kpc)
rho01 = 1.64e7          # core density (Msun/kpc**3)
rtrunc1 = 300.          # truncation (virial) radius
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)       # potential (km/s)^2
x1_0 = 0.               # positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.              # velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sgr
Mass2_scale = 2e10
Potential2 = funcs.Kepler_Potential(Mass2_scale)       # potential (km/s)^2
x2_0 = -20.             # positions (kpc)
y2_0 = 0.
z2_0 = 0
vx2_0 = 0.            # velocities (km/s):
vy2_0 = 0.
vz2_0 = -330.
