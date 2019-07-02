'''
Author: Thorsten Tepper-Garcia
Date: 02/07/2019
'''

import config.phys_consts as pc
from utils import funcs

# Parameters appropriate for the Solar circle in the Galaxy
Mass = 1.0e11													# mass (Msun)
Potential = funcs.Kepler_Potential(amp=pc.Grav*Mass)			# potential (km/s)^2
t_0 = 0.														# initial time (Gyr)
t_1 = 10.														# total time (in time units)
x_0 = 8.														# positions (kpc)
y_0 = 0.
vx_0 = 0.														# velocities (km/s):
vy_0 = funcs.v_circ(x_0,y_0, amp = pc.Grav*Mass)				# circular velocity at r_0
