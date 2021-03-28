'''
Author: Thorsten Tepper Garcia
Date: 29/04/2020

Orbit intended for MW-Sgr interaction with AGAMA paper


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0             # initial time (Gyr)
t_1 = 3.0e-1           # total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-4        # integration time step

# Milky Way
rs1 = 1.5e1             # NFW scale radius (kpc)
rho01 = 1.57e7          # core density (Msun/kpc**3)
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
rs2 = 1.0e1             # NFW scale radius (kpc) -> CORRECT?
rho02 = 8.6e6            # core density (Msun/kpc**3) -> CORRECT?
rtrunc2 = 25.          # truncation (virial) radius -> CORRECT?
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.NFW_Potential(Mass2_scale,rs2)       # potential (km/s)^2
x2_0, y2_0, z2_0 = 14.83024,0.00000,50.09880
vx2_0, vy2_0, vz2_0 = -134.24384,-0.00000,-130.87080

# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=3.e-1*rs2)	# DF function
