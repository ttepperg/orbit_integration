'''
Author:	Thorsten Tepper Garcia
Date:	27/06/2019

This is just a test run.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e1														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 1.2e1														# NFW concentration
rs1 = 2.0e1														# NFW scale radius (kpc)
rho01 = 9.0e6													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)						# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)										# total ('virial') mass (Msun)
x1_0 = 7.1286E-09													# positions (kpc)
y1_0 = 1.2991E-08
z1_0 = -2.4789E-08
vx1_0 = -1.2330E-09												# velocities (km/s):
vy1_0 = -4.2230E-09
vz1_0 = 2.7000E-09

# Large Magellanic Clouds
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2)							# 'cumulative' mass, trivially equal to Mass2
x2_0 = -2.8085E+01												# positions (kpc)
y2_0 = 3.3417E+01
z2_0 = 1.6563E+02
vx2_0 = 1.7074E+01												# velocities (km/s):
vy2_0 = 8.6930E+01
vz2_0 =  -1.4530E+01


#output
print("Virial radius of body 1: {}".format(Rvir1))
