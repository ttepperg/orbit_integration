'''
Author:	Thorsten Tepper Garcia

This setting is essentially identical to lmc_mw_infall_galpy.py, but takes
dynamical friction by the MW onto the LMC into account. See lmc_mw_infall_galpy.py for
more information.

The goal is to roughly match the LMC's orbit given by Galpy found in:

https://galpy.readthedocs.io/en/v1.4.0/orbit.html?highlight=lmc#new-in-v1-4-example-the-orbit-of-the-large-magellanic-cloud-in-the-presence-of-dynamical-friction

with the softening length of the LMC is set to 5 kpc.

In the Galpy example, the LMC by a point-like object, as is done here. A notable difference
between Galpy's example and this calculation however, is that Galpy assumes the MW to be fixed
at all times, while here both the MW and the LMC move in response to each other's gravitational
pull (although the MW's motion is small given the large MW : LMC mass ratio). 

In addition, Galpy uses a somewhat different parametrization for the dynamical friction.

Thus, some small differences between the two results may be expected.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 1.53e1														# NFW concentration
rs1 = 1.6e1														# NFW scale radius (kpc)
rho01 = 1.307e7													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)						# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)										# total ('virial') mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
Mass2 = 5.0e10													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
x2_0 = 0.459													# positions (kpc)
y2_0 = -41.288
z2_0 = -27.149
vx2_0 = 39.708													# velocities (km/s):
vy2_0 = -241.373
vz2_0 = 230.308


# Dynamical friction settings
soft_length2 = 5.0													# softening length of LMC
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Info
# print("Virial radius (kpc): {:.3f}".format(Rvir1))

