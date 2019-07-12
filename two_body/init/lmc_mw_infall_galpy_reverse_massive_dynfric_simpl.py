'''
Author:	Thorsten Tepper Garcia

This setup is in essence identical to lmc_mw_infall_galpy_reverse_massive.py, but
including dynamical friction adopting a simplified formula.


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.02e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 1.4e1														# NFW concentration
rs1 = 1.6e1														# NFW scale radius (kpc)
rho01 = 1.26e7													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)						# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)										# total ('virial') mass (Msun)
x1_0 = 8.1707E-09												# positions (kpc)
y1_0 = 6.7203E-09
z1_0 = -3.4977E-08
vx1_0 = -7.5039E-10												# velocities (km/s):
vy1_0 = -3.8106E-09
vz1_0 = 6.4645E-10
Dynamical_Friction1 = funcs.dyn_friction_simpl()				# dynamical friction function

# Large Magellanic Clouds
Mass2 = 1.0e10													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
x2_0 = -2.1025E+01												# positions (kpc)
y2_0 = -4.6087E+01
z2_0 = 6.6870E+01
vx2_0 = 4.7124E+01												# velocities (km/s):
vy2_0 = 2.4655E+02
vz2_0 = -3.4774E+01

