'''
Author:	Thorsten Tepper Garcia

This set of parameters is intended to reproduce the infall of the Large Magellanic Cloud (LMC)
onto the Milky Way (MW), using backwards integration.

The LMC position and velocity are taken from Galpy
(see https://galpy.readthedocs.io/en/v1.4.0/orbit.html?highlight=lmc#new-in-v1-4-example-the-orbit-of-the-large-magellanic-cloud-in-the-presence-of-dynamical-friction):


LMC initial orbital parameters:
Galactocentric position (x,y,z) = (0.459,-41.288,-27.149) kpc
Galactocentric velocity (vx,vy,vz) = (39.708,-241.373,230.308) km/s

Note that the sing of the velocity components are different from other sources,e.g.
https://arxiv.org/pdf/1301.0832.pdf, likely to the difference in the sign convention
for the proper motion along the right ascension. Galpy retrieves the values
directly from the SIMBAD database.

The MW is approximated by a static NFW DM halo; the LMC by a point-like object.
The initial MW DM halo parameters were taken from galpy's MWPotential2014 (see Bovy (2014; their
table 1, at https://arxiv.org/pdf/1412.3451.pdf), but the MW DM halo's centra density has been
increased (since the other components are ignored) so as to roughly match the LMC's orbit given
by Galpy found in:

https://galpy.readthedocs.io/en/v1.4.0/orbit.html?highlight=lmc#new-in-v1-4-example-the-orbit-of-the-large-magellanic-cloud-in-the-presence-of-dynamical-friction


Note that this model does not take the presence of the Small Magellanic Cloud (SMC)
or dynamical friction into account (yet).


'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = -1.0226e1													# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Milky Way
c1 = 1.53e1														# NFW concentration
rs1 = 1.6e1														# NFW scale radius (kpc)
rho01 = 1.392e7													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)						# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)										# total ('virial') mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Large Magellanic Clouds
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
x2_0 = 0.459														# positions (kpc)
y2_0 = -41.288
z2_0 = -27.149
vx2_0 = 39.708													# velocities (km/s):
vy2_0 = -241.373
vz2_0 = 230.308


# Info
# print("Virial radius (kpc): {}".format(Rvir1))
