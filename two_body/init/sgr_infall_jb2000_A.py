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


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.135e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
rho01 = 2.e8
a1 = 1.e0														# PITS scale radius
Potential1 = funcs.PITS_Potential(rho0=rho01,a=a1)				# potential (km/s)^2
Mass1_cum = funcs.PITS_Mass(rho0=rho01,a=a1)					# mass function
rt1 = [0.,0.,1000.]												# truncation radius
Mass1 = Mass1_cum(*rt1)											# total mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
rho02 = 1.2e8
a2 = 1.e0														# PITS scale radius
Potential2 = funcs.PITS_Potential(rho0=rho02,a=a2)				# potential (km/s)^2
Mass2_cum = funcs.PITS_Mass(rho0=rho02,a=a2)					# mass function
rt2 = [0.,0.,70.]												# truncation radius
Mass2 = Mass2_cum(*rt2)											# total mass (Msun)
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
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function

# Info
# print("Mass of MW [Msun]: {:E}".format(Mass1))
# print("Mass of Sgr [Msun]: {:E}".format(Mass2))
# 
# # Check
# import math
# sigma2 = funcs.PITS_VelDisp(rho01,a1)
# print(math.sqrt(sigma2(*rt2)))
# print(funcs.PITS_Vinf(rho01,a1))
# exit()
