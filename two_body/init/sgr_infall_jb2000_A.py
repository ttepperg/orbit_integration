'''
Author:	Thorsten Tepper Garcia

This is setting is an attempt to reproduce model A of Jiang & Binney (2000), with
the following differences:
- MW mass profile (see their equation 2); which affects dynamical friction
- Sgr mass profile
- Absence of Galactic disc (see their equation 1a,b)
- Absence of mass loss experienced by Sgr due to tidal stripping (see their equation 4a,b)
- Differences in the calculation of the dynamical friction:
	- Variable Coloumb logarithm  (they adopt a constant value Log L = 8.5)

Therefore, some differences between theirs and our results are expected.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.135e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1 = 2.5e12													# total mass (Msun)
rs1 = 3.e1														# Plummer scale radius
Potential1 = funcs.Plummer_Potential(mass=Mass1,a=rs1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
Mass2 = 1.0e10													# total mass (Msun)
rs2 = 1.e-2														# Plummer scale radius
Potential2 = funcs.Plummer_Potential(mass=Mass2,a=rs2)			# potential (km/s)^2
x2_0 = 0.														# positions (kpc)
y2_0 = 0.
z2_0 = -250.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = -103.
vz2_0 = 0.


# Dynamical friction settings
soft_length2 = 1.0													# softening length of Sgr
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Info
# print("Virial radius (kpc): {}".format(Rvir1))
