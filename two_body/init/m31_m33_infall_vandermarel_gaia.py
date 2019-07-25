'''
Author:	Thorsten Tepper Garcia

This set of parameters is intended to reproduce the semi-analytic calculation
for the infall of M33 onto M31 using the most recent Gaia proper motion measurements
for both galaxies by van der Marel et al. (2019), available from:
https://ui.adsabs.harvard.edu/abs/2019ApJ...872...24V/abstract

To have a direct comparison with previous work (Semczuck et al. 2018), we follow
these assumptions:

M31 is approximated by a NFW DM halo with a virial mass of 2x10^12 Msun.
M33 is approximated by a NFW DM halo with a virial mass of 5x10^11 Msun.

The Galactocentric velocity of M31 from the Gaia DR2 data alone is:

v_M31 = (0 ± 75, −176 ± 51, −84 ± 73) km/s

The Galactocentric velocity of M33 from the Gaia DR2 data alone is:

v_M33 = (49 ± 74, 14 ± 70, 28 ± 73) km/s

For now, we neglect the error bars in both the relative coordinates and velocities
and adopt the central measured values only.

The velocity of M33 relative to M31 is thus

v_rel = (49, 190, 112)

The orbit is integrated for 10 Gyr to allow for a broader history.

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.022e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 11.65														# NFW scale radius (kpc)
rho01 = 4.2e7													# core density (Msun/kpc**3)
Mass1 = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1,rs1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
rs2 = 17.88														# NFW scale radius (kpc)
# rho02 = 3.88e6													# core density (Msun/kpc**3)
rho02 = 4.45e6													# core density (Msun/kpc**3)
Mass2 = 4. * pi * rho02 * rs2**3
Potential2 = funcs.NFW_Potential(Mass2,rs2)			# potential (km/s)^2
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = 49.													# velocities (km/s):
vy2_0 = 190.
vz2_0 = 112.


# Dynamical friction settings
soft_length2 = 28.5													# softening length of M33
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=soft_length2)	# dynamical friction function

