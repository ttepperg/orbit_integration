'''
Author:	Thorsten Tepper Garcia

This setting essentially corresponds to lmc_mw_infall_galpy_dynfric.py,
but integrated forward in time, using the state vectors provided by
the backwards integration as initial conditions.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 1.53e1														# NFW concentration
rs1 = 1.6e1														# NFW scale radius (kpc)
rho01 = 1.307e7												# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)			# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)									# total ('virial') mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
Mass2 = 5.0e10													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)		# potential (km/s)^2
x2_0 = -44.484													# positions (kpc)
y2_0 = 13.772
z2_0 = -456.733
vx2_0 = 12.913													# velocities (km/s):
vy2_0 = -52.528
vz2_0 = 95.000

# Dynamical friction settings
soft_length2 = 5.0																	# softening length of LMC
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Info
# print("Virial radius (kpc): {:.3f}".format(Rvir1))

