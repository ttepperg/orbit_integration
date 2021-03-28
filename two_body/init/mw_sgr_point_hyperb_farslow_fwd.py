'''
Author: Thorsten Tepper Garcia
Date: 29/04/2020

Modifications
17/05/2020 - changed satellite mass from 5.6e10 to 2e10
			- changed name from mw_sgr_hyperb_fwd_point_farslow.py to
				mw_sgr_point_hyperb_farslow_fwd.py
			- swapped y2_0 and z2_0 coordinates
			- changed rho01 from 1.57e7 to 1.64e7

THIS is the orbit actually used to model the MW-Sgr interaction with AGAMA paper (Bland-Hawthorn & Tepper-Garcia 2021)

Note that t_1 here is larger than in mw_sgr_point_hyperb_farslow_bck.py to verify whether the orbit is truly hyperbolic (appers to be not, which is why in our paper we introduced a mass-reduction scheme). To check that both the forward and backward integrations are consistent with one another. set t_1 = 0.1 here and compare the present-day orbital parameters with the infall parameters in mw_sgr_point_hyperb_farslow_bck.py.


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0             # initial time (Gyr)
t_1 = 6.0e-1           # total time (time unit ~ 0.978 Gyr)
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
x2_0, y2_0, z2_0 = -11.09052,-0.00000,28.59071
vx2_0, vy2_0, vz2_0 = -145.41697,0.00000,-220.22636


#EOF
