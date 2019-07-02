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
semimajor_a = 2.												# must be > r_0/2
vy_0 = funcs.vis_viva(x_0,y_0,amp=pc.Grav*Mass,a=semimajor_a)	# elliptical orbit with semi-major axis a > r/2
