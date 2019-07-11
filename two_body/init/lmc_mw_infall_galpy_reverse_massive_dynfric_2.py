'''
Author:	Thorsten Tepper Garcia
Date:	10/07/2019

This setup is in essence identical to lmc_mw_infall_galpy_reverse_massive_dynfrict.py, but
with the body indices swapped to check the invariance of the code with respect to
the order of the initial conditions.


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.02e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c2 = 1.4e1														# NFW concentration
rs2 = 1.6e1														# NFW scale radius (kpc)
rho02 = 1.26e7													# core density (Msun/kpc**3)
Rvir2 = c2*rs2
Potential2 = funcs.NFW_Potential(rho02,rs2)						# potential (km/s)^2
Mass2_cum = funcs.NFW_Mass(rho02,rs2)
Mass2 = Mass2_cum(Rvir2)										# total ('virial') mass (Msun)
x2_0 = 8.1707E-09												# positions (kpc)
y2_0 = 6.7203E-09
z2_0 = -3.4977E-08
vx2_0 = -7.5039E-10												# velocities (km/s):
vy2_0 = -3.8106E-09
vz2_0 = 6.4645E-10
Dynamical_Friction2 = funcs.dyn_friction_simpl					# dynamical friction function

# Large Magellanic Clouds
Mass1 = 1.0e10													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(amp=pc.Grav*Mass1)			# potential (km/s)^2
x1_0 = -2.1025E+01												# positions (kpc)
y1_0 = -4.6087E+01
z1_0 = 6.6870E+01
vx1_0 = 4.7124E+01												# velocities (km/s):
vy1_0 = 2.4655E+02
vz1_0 = -3.4774E+01

