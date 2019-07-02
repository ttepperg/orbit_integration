import config.phys_consts as pc
from utils import funcs

t_0 = 0.														# initial time (Gyr)
t_1 = 10.														# total time (time unit ~ 0.978 Gyr)

# Body 1
Mass1 = 1.0e6													# mass (Msun)
Potential1 = funcs.Kepler_Potential(amp=pc.Grav*Mass1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.

# Body 2
Mass2 = 1.0e6													# mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 2.
