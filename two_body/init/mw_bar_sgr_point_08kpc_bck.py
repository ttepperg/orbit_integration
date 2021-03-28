'''
AUTHOR:
	Thorsten Tepper Garcia

CREATED: 29/03/2021

ABOUT:
	This is the FAR orbit intended to model the MW(bar)-Sgr interaction with AGAMA. It is based on model hbd_9

NOTES:
	- The mass model of the host galaxy must be consistent with its *total* mass (i.e. taking all components into account). The best way to achieve this is first to try adjust the parameters of the (dominant) DM halo (rs1, rho01,rtrunc1) to match the total mass within the virial (truncation) radius. Then adjust the density normalization to match the actual (total) mass enclosed within the orbit of the satellite, i.e. setting rtrunc1 =  radius of orbit at crossing point. Generally, this means to slightly increase rho01 beyond its value for a single-component DM halo. Note that this will lead to a higher total mass out to the virial radius, which is fine.

	For model hbd_9, the enclosed mass within 20 kpc is ~2x10^11 Msun. The total mass out to 300 kpc is ~1.23x10^12 Msun.

	DO NOT FORGET to set 'rtrunc1 = virial radius' of the host galaxy!


CODE CALL & OUTPUT:

	mac-ina:two_body tepper$ py37 two_body_orbit.py mw_bar_sgr_20kpc_bck.py

	Constants and units:

	  Gravitational constant:   4.301E-06 kpc km^2 / Msun s^2
	               Mass unit:       1.000 solMass = 1.988410E+33 g
	             Length unit:       1.000 kpc = 3.085678E+21 cm
	           Velocity unit:       1.000 km / s = 1.000000E+05 cm / s
	               Time unit:       0.978 Gyr = 3.085678E+16 s

	Removed unnecessary extension .py in input parameter file.

	Gathering input parameters from file mw_bar_sgr_20kpc_bck...
	Done.

	WARNING:
		orbital circumference not defined for e >= 1.


	Initial (osculating) orbital parameters of the system:

	                    Potential of body 1:        NFW_Pot
	                    Potential of body 2:     Kepler_Pot
	                   Bound mass of body 1:     1.7828E+12
	                   Bound mass of body 2:     2.0000E+10

	          Initial rel. distance (r21_0):        20.0000
	             Initial rel. speed (v21_0):       330.0000
	      Initial tangential vel. (v_tan_0):       330.0000
	          Initial radial vel. (v_rad_0):         0.0000
	   Rel. specific ang. mom. vec. (h_vec): (-0.000,-6600.000,-0.000)
	                           [normalised]: (-0.000,-1.000,-0.000)
	                         [along z-axis]: (0.000,-0.000,6600.000)
	            Rel. specific ang. mom. (h):      6600.0000

	                       Eccentricity (e):         1.3488
	                  Semi-latus rectum (p):        46.9764
	                     Semimajor axis (a):       -57.3363
	                     Semiminor axis (b):        51.8985
	       Orbital inclination (psi_0; deg):        90.0000
	                 Ascending node (n_vec): (6600.000,0.000,0.000)
	      Long. of asc. node (Omega_0; deg):         0.0000
	            Eccentricity vector (e_vec): (-1.349,0.000,0.000)
	             Apsidal angle (phi_0; deg):       180.0000
	   Argument of periapsis (omega_0; deg):       180.0000
	                   Rel. pericentre (rp):        20.0000
	                Vel. at pericentre (vp):       330.0000
	     Rel. specific potential energy (V):    -1.5092E+05
	       Rel. specific kinetic energy (T):     5.4450E+04
	         Rel. specific total energy (E):    -9.6472E+04


	Orbit integration (backwards):
		Time range [t0,t1] = [-0.1,0.0]
		Time steps 1000;  Step size: 1.00E-04

	WARNING:
		orbital circumference not defined for e >= 1.

		Maximum time step to avoid energy drift:          NAN
		Expected accumulated error in x and v of order 1.00E-10

	WARNING:
		orbital circumference not defined for e >= 1.


	Output with timestep frequency output_freq=10 written to file:

		./output/mw_bar_sgr_20kpc_bck_out.dat

	Conservation laws:
	                          Energy conservation to better than 8.498E-06 %.
	    Angular momentum conservation (magnitude) to better than 1.589E-10 %.
	    Angular momentum conservation (direction) to better than 0.000E+00 %.


	Initial / final state vectors:
		Infall relative position: (-10.87381,0.00000,28.44463)
		Infall relative distance: 30.45221
		Infall relative velocity: (-149.58385,-0.00000,-215.66883)
		Infall relative speed: 262.46595
		Present-day relative position: (-20.00000,0.00000,0.00000)
		Present-day relative distance: 20.00000
		Present-day relative velocity: (0.00000,0.00000,-330.00000)
		Present-day relative speed: 330.00000

	The following are only relevant if masses are allowed to evolved:
		Infall bound mass of body 1: 1.782763E+12
		Infall bound mass of body 2: 2.000000E+10

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0             # initial time (Gyr)
t_1 = -1.0e-1           # total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-4        # integration time step

# Milky Way
rs1 = 1.9e1             # NFW scale radius (kpc)
rho01 = 1.1e7          # core density (Msun/kpc**3)
rtrunc1 = 300.          # truncation (virial) radius
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)       # potential (km/s)^2
x1_0 = 0.               # positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.              # velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Sgr
Mass2_scale = 2e10
Potential2 = funcs.Kepler_Potential(Mass2_scale)       # potential (km/s)^2
x2_0 = -20.             # positions (kpc)
y2_0 = 0.
z2_0 = 0
vx2_0 = 0.            # velocities (km/s):
vy2_0 = 0.
vz2_0 = -330.
