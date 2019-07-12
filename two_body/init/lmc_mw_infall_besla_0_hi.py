'''
Author:	Thorsten Tepper Garcia

This set of parameters is intended to reproduce the infall of the Large Magellanic Cloud (LMC)
onto the Milky Way (MW), using backwards integration, i.e. inverting the sign of the
observed velocity.

The LMC parameters correspond to the K1 mean values as given by Besla et al. (2007, their
table 1).

The MW is approximated by a static NFW halo; the LMC by a point-like object.

The MW DM halo parameters are taken from Besla et al. (2010, their table 3).
This model corresponds to the 'high-mass' MW model.

Note that this model does not take the presence of the Small Magellanic Cloud (SMC)
or dynamical friction into account (yet).


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e1														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 9.0e0														# NFW concentration
rs1 = 3.59e1													# NFW scale radius (kpc)
rho01 = 2.5e6													# core density (Msun/kpc**3)
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
Mass2 = 2.0e10													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2)							# 'cumulative' mass, trivially equal to Mass2
x2_0 = -0.8														# positions (kpc)
y2_0 = -41.5
z2_0 = -26.9
vx2_0 = -86.													# velocities (km/s):
vy2_0 = -268.
vz2_0 = 252.


#output
print("Virial radius of Milky Way: {}".format(Rvir1))
