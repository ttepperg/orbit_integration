'''
Author: Thorsten Tepper-Garcia
Date: 02/07/2019
'''

import config.phys_consts as pc
from utils import funcs

Mass = 1.0e6													# mass (Msun)
Potential = funcs.Kepler_Potential(amp=pc.Grav*Mass)			# potential (km/s)^2
t_0 = 0.														# initial time (Gyr)
t_1 = 10.														# total time (in time units)
x_0 = 1.														# positions (kpc)
y_0 = 0.
vx_0 = 0.														# velocities (km/s):
vy_0 = 1.
