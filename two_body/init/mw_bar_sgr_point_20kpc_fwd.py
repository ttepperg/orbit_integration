'''
AUTHOR:
	Thorsten Tepper Garcia

CREATED: 29/03/2021

ABOUT:
	This is the FAR orbit intended to model the MW(bar)-Sgr interaction with AGAMA. It is based on model hbd_9

NOTES:
	The initial orbital parameters are obtained from the infall parameters calculated with mw_bar_sgr_point_20kpc_bck.py. Recall that a self-consistent result implies infall parameters calcualted with one init file (e.g. bck) matching the present-day parameters of its counterpart (e.g. fwd).

	Set t_1 = 6.0e-1 to integrate beyond the crossing point at z=0 and R=20 to check whether the orbit is truly hyperbolic (appers to be not).

CODE CALL & OUTPUT:

	mac-ina:two_body tepper$ py37 two_body_orbit.py mw_bar_sgr_point_20kpc_fwd.py

	Constants and units:

	  Gravitational constant:   4.301E-06 kpc km^2 / Msun s^2
	               Mass unit:       1.000 solMass = 1.988410E+33 g
	             Length unit:       1.000 kpc = 3.085678E+21 cm
	           Velocity unit:       1.000 km / s = 1.000000E+05 cm / s
	               Time unit:       0.978 Gyr = 3.085678E+16 s

	Removed unnecessary extension .py in input parameter file.

	Gathering input parameters from file mw_bar_sgr_point_20kpc_fwd...
	Done.

	Initial (osculating) orbital parameters of the system:

	                    Potential of body 1:        NFW_Pot
	                    Potential of body 2:     Kepler_Pot
	                   Bound mass of body 1:     1.7828E+12
	                   Bound mass of body 2:     2.0000E+10

	          Initial rel. distance (r21_0):        30.4522
	             Initial rel. speed (v21_0):       262.4659
	      Initial tangential vel. (v_tan_0):       216.7331
	          Initial radial vel. (v_rad_0):      -148.0377
	   Rel. specific ang. mom. vec. (h_vec): (0.000,-6599.999,0.000)
	                           [normalised]: (0.000,-1.000,0.000)
	                         [along z-axis]: (0.000,-0.000,6599.999)
	            Rel. specific ang. mom. (h):      6599.9991

	                       Eccentricity (e):         0.6628
	                  Semi-latus rectum (p):        29.5195
	                     Semimajor axis (a):        52.6516
	                     Semiminor axis (b):        39.4240
	       Orbital inclination (psi_0; deg):        90.0000
	                 Ascending node (n_vec): (6599.999,0.000,-0.000)
	      Long. of asc. node (Omega_0; deg):         0.0000
	            Eccentricity vector (e_vec): (-0.608,0.000,-0.265)
	             Apsidal angle (phi_0; deg):      -156.4307
	   Argument of periapsis (omega_0; deg):       203.5693
	                   Rel. pericentre (rp):        17.7526
	                Vel. at pericentre (vp):       371.7769
	                    Rel. apocentre (ra):        87.5506
	                 Vel. at apocentre (va):        75.3850
	                Rel. orbital period (T):         1.9761
	        Approx. orbit circumference (u):       290.7563
	        Approx. pericentric period (Tp):         0.7821
	     Rel. specific potential energy (V):    -1.3092E+05
	       Rel. specific kinetic energy (T):     3.4444E+04
	         Rel. specific total energy (E):    -9.6472E+04


	Orbit integration (forward):
		Time range [t0,t1] = [0.0,0.1]
		Time steps 1000;  Step size: 1.00E-04
		Maximum time step to avoid energy drift: 1.760281E-01
		Expected accumulated error in x and v of order 1.00E-10

	Output with timestep frequency output_freq=10 written to file:

		./output/mw_bar_sgr_point_20kpc_fwd_out.dat

	Conservation laws:
	                          Energy conservation to better than 8.498E-06 %.
	    Angular momentum conservation (magnitude) to better than 2.088E-10 %.
	    Angular momentum conservation (direction) to better than 0.000E+00 %.


	Initial / final state vectors:
		Infall relative position: (-10.87381,0.00000,28.44463)
		Infall relative distance: 30.45220
		Infall relative velocity: (-149.58385,-0.00000,-215.66883)
		Infall relative speed: 262.46595
		Present-day relative position: (-20.00000,0.00000,-0.00000)
		Present-day relative distance: 20.00000
		Present-day relative velocity: (0.00002,0.00000,-330.00000)
		Present-day relative speed: 330.00000

	The following are only relevant if masses are allowed to evolved:
		Present-day bound mass of body 1: 1.782763E+12
		Present-day bound mass of body 2: 2.000000E+10

'''

import config.phys_consts as pc
from utils import funcs
from math import pi


t_0 = 0.0e0             # initial time (Gyr)
t_1 = 1.0e-1           # total time (time unit ~ 0.978 Gyr)
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
x2_0, y2_0, z2_0 = -10.87381,0.00000,28.44463
vx2_0, vy2_0, vz2_0 = -149.58385,-0.00000,-215.66883


#EOF
