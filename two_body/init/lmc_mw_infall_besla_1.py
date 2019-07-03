'''
Author:	Thorsten Tepper Garcia
Date:	03/07/2019

Experimental

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e1														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 1.2e1														# NFW concentration
rs1 = 2.0e1														# NFW scale radius (kpc)
rho01 = 7.0e6													# core density (Msun/kpc**3)
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
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2)							# 'cumulative' mass, trivially equal to Mass2
x2_0 = -1.														# positions (kpc)
y2_0 = -41.
z2_0 = -28.
vx2_0 = 79.														# velocities (km/s):
vy2_0 = 227.
vz2_0 = -208.


#output
print("Virial radius of body 1: {}".format(Rvir1))
