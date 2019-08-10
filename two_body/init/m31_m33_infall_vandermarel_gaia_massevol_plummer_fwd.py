'''
Author:	Thorsten Tepper Garcia

This setup is essentially aims at tuning the infall of M33 onto M31
iteratively by hand, alternating between this setup and
m31_m33_infall_vandermarel_gaia_massevol.py.

The boundary conditions are:
- present-day position and velocity of M33 relative to M31
- M33's present-day mass => ~3x10^11 Msun (at the low end of the range
  estimated by Corbelli et al. 2014)

Initial conditions: the most recent Gaia proper motion measurements
for both galaxies by van der Marel et al. (2019), available from:
https://ui.adsabs.harvard.edu/abs/2019ApJ...872...24V/abstract

The Galactocentric velocity of M31 from the Gaia DR2 data alone is:

v_M31 = (0 ± 75, −176 ± 51, −84 ± 73) km/s

The Galactocentric velocity of M33 from the Gaia DR2 data alone is:

v_M33 = (49 ± 74, 14 ± 70, 28 ± 73) km/s

For now, we neglect the error bars in both the relative coordinates and velocities
and adopt the central measured values only.

The velocity of M33 relative to M31 is thus

v_rel = (49, 190, 112)

The orbit is integrated for 10 Gyr to allow for a broader history.

The mass model for M31 corresponds to the high-mass end adopted by Patel et al. (2017).
The mass model for M33 is based on Corbelli et al. (2014).


The present-day state vectors are:

	Present-day relative position: (-95.74430,-122.02125,-128.66744)
	Present-day relative distance: 201.52287
	Present-day relative velocity: (45.85365,188.81194,108.79276)
	Present-day relative speed: 222.68446


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.33e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 11.65														# NFW scale radius (kpc)
rho01 = 4.2e7													# core density (Msun/kpc**3)
rtrunc1 = 326.													# truncation (virial) radius
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)				# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
rs2 = 17.88														# NFW scale radius (kpc)
rho02 = 6.135e6												# core density (Msun/kpc**3)
rtrunc2 = 197.													# truncation (virial) radius
Mass2_scale = 4. * pi * rho02 * rs2**3
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)				# potential (km/s)^2
x2_0 = -944.23627												# positions (kpc)
y2_0 = 444.40255
z2_0 = -672.72794
vx2_0 = 81.41245												# velocities (km/s):
vy2_0 = -13.60903
vz2_0 = 66.94253


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)		# dynamical friction function

# Mass evolution (stripping)
Mass1_cum = funcs.NFW_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Plummer_Mass(Mass2_scale,rs2)
Mass2_evol = funcs.mass_bound(m1_func=Mass1_cum,m2_func=Mass2_cum)	# mass evolution function
