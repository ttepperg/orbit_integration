'''
Author:	Thorsten Tepper Garcia

This set of parameters is intended to reproduce the semi-analytic calculation
for the infall of M33 onto M31 using the most recent Gaia proper motion measurements
for both galaxies by van der Marel et al. (2019), available from:
https://ui.adsabs.harvard.edu/abs/2019ApJ...872...24V/abstract

with the following differences:

- We approximate M31 by a DM halo only, and ignore additional
component (e.g. stellar bulge, stellar disc). The mass of these components
is added in each case to the host DM halo.
- The expression for the Couloumb logarithm differs slightly, more specifically,
the constant factor in the denominator we use is 1.4 rather than 1.2


Initial conditions:

This set of parameters is intended to reproduce the semi-analytic calculation
for the infall of M33 onto M31 using the most recent Gaia proper motion measurements
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

The mass models for the galaxies are based on the extreme cases explored
in Patel et al. (2017). Here, we adopt:

M31 -> high mass case (model b in that paper)
M33 -> high mass case


'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.33e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# M31
rs1 = 35.15														# NFW scale radius (kpc)
rho01 = 2.56e6													# core density (Msun/kpc**3)
Mass1_scale = 4. * pi * rho01 * rs1**3
rtrunc1 = 329.													# truncation (virial) radius
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)				# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# M33
rs2 = 21.														# scale radius (kpc)
Mass2_scale = 2.5e11													# mass scale (Msun)
Potential2 = funcs.Plummer_Potential(Mass2_scale,rs2)		# potential (km/s)^2
x2_0 = -97.2													# positions (kpc)
y2_0 = -121.6
z2_0 = -129.8
vx2_0 = 49.													# velocities (km/s):
vy2_0 = 190.
vz2_0 = 112.


# Dynamical friction settings
Dynamical_Friction1 = funcs.dyn_friction_maxwell(eps=rs2)	# dynamical friction function
