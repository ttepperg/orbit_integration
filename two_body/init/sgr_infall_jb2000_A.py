'''
Author:	Thorsten Tepper Garcia

This setup is an attempt to reproduce model A of Jiang & Binney (2000), with
the following differences:

- Both the MW and Sgr are free to move
- Our MW and Sgr mass profiles do not have a tapper function (see their equation 2), which
  slightly changes the potential
- Absence of dilution factor in mass loss experienced by Sgr due to tidal stripping (see their equation 4a,b)
- Differences in the calculation of the dynamical friction:
	- Variable Coloumb logarithm  (they adopt a constant value Log L = 8.5)

Therefore, some differences between theirs and our results are expected.

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.135e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
a1 = 1.e0														# PITS scale radius
rho01 = 2.e8
Mass1_scale_scale = 4. * pi * rho01 * a1**3
Potential1 = funcs.PITS_Potential(Mass1_scale_scale,a1)						# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
a2 = 1.e0														# PITS scale radius
rho02 = 1.2e8
Mass2_scale_scale = 4. * pi * rho02 * a2**3
Potential2 = funcs.PITS_Potential(Mass2_scale_scale,a2)						# potential (km/s)^2
x2_0 = 0.														# positions (kpc)
y2_0 = 0.
z2_0 = -250.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = -103.
vz2_0 = 0.


# Dynamical friction settings
soft_length2 = 2.0e1												# softening length of Sgr
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

# Mass loss
Mass1_cum = funcs.PITS_Mass(Mass1_scale_scale,a1)							# mass function
Mass2_cum = funcs.PITS_Mass(Mass2_scale_scale,a2)							# mass function
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function

# Info
rt1 = [0.,0.,1000.]												# truncation radius
m1 = Mass1_cum(*rt1)											# total mass (Msun)
rt2 = [0.,0.,70.]												# truncation radius
m2 = Mass2_cum(*rt2)											# total mass (Msun)
print("Mass of MW [Msun]: {:E}".format(m1))
print("Mass of Sgr [Msun]: {:E}".format(m2))
r0 = [x2_0,y2_0,z2_0]
rt = funcs.tidal_radius(*r0,m1_func=Mass1_cum,m2_func=Mass2_cum)
print(rt)
print("{:E}".format(Mass2_cum(rt)))
# exit()
# 
# # Check
# import math
# sigma2 = funcs.PITS_VelDisp(rho01,a1)
# print(math.sqrt(sigma2(*rt2)))
# print(funcs.PITS_Vinf(rho01,a1))
# exit()
