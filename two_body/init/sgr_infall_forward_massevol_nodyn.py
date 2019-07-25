'''
Author:	Thorsten Tepper Garcia

This setup is essentially identical to sgr_infall_forward_massevol.py, but
neglecting dynamical friction.

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
Mass1 = 1.325e12												# total mass (Msun)
a1 = 38.35														# scale radius (kpc)
Potential1 = funcs.Hernquist_Potential(mass=Mass1,a=a1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sagittarius dwarf
Mass2 = 1.3e10													# total mass (Msun)
a2 = 9.81														# scale radius (kpc)
Potential2 = funcs.Hernquist_Potential(mass=Mass2,a=a2)			# potential (km/s)^2
x2_0 = -259.36180422351384									# positions (kpc)
y2_0 = -40.32912811156189
z2_0 = 125.16118908506236
vx2_0 = 65.79976992535424									# velocities (km/s):
vy2_0 =  13.89047520759926
vz2_0 = -71.08450826725453

# Dynamical friction settings
# soft_length2 = 1.0																	# softening length of LMC
# Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function


# Mass evolution
Mass1_cum = funcs.Hernquist_Mass(mass=Mass1,a=a1)
Mass2_cum = funcs.Hernquist_Mass(mass=Mass2,a=a2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function

