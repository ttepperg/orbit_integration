'''
Author:	Thorsten Tepper Garcia
Date:	11/07/2019

Experimental!

This is essentially identical to lmc_mw_infall_plummer.py but includes
the effect of dynamical friction of the MW onto the LMC.


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.02e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1 = 8.e11
rs = 1.6e1														# NFW scale radius (kpc)
Potential1 = funcs.Plummer_Potential(amp=pc.Grav*Mass1,a=rs)	# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.
Dynamical_Friction1 = funcs.dyn_friction_maxwell(pot=Potential1,eps=1.)	# dynamical friction function

# Large Magellanic Clouds
Mass2 = 1.0e10													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = -1.														# positions (kpc)
y2_0 = -41.
z2_0 = -28.
vx2_0 = 79.														# velocities (km/s):
vy2_0 = 227.
vz2_0 = -208.

# test area
# import math
# r0 = [30.,0.5,1.]
# sigma2 = funcs.Plummer_VelDisp(amp=pc.Grav*Mass1,a=rs)
# print(math.sqrt(sigma2(*r0)))
# 
# exit()
