'''
Author: Thorsten Tepper Garcia
Date: 23/05/2020

Copy of mw_sgr_point_hyperb_farslow_df_fwd.py but with stronger DF

Orbit intended for MW-Sgr interaction with AGAMA paper


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0             # initial time (Gyr)
t_1 = 6.0e-1            # total time (time unit ~ 0.978 Gyr)
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
Mass2_scale = 2.e10
Potential2 = funcs.Kepler_Potential(Mass2_scale)       # potential (km/s)^2
x2_0, y2_0, z2_0 = -10.91004,-0.00000,29.98046
vx2_0, vy2_0, vz2_0 = -149.21618,0.00000,-245.12100


# Dynamical friction settings
# Note: the value of eps has been found to be a good match to the
# orbit obtained from an Nbody model ran using the orbital inital
# parameters given in mw_sgr_point_hyperb_farslow_df_fwd and comparing
# to the orbit resulting from mw_sgr_point_hyperb_farslow_fwd
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=2.7e-1)	# DF function


#EOF
