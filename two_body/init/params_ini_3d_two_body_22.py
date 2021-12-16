'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved equal mass, two-body system, with a point mass and an extended body on an elliptic orbit with the centre of attraction at its centre rather than at one focus, characterised by the following parameters:

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.0000
      Initial tangential vel. (v_tan_0):         1.0000
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,1.000)
                           [normalised]: (0.000,0.000,1.000)
                         [along z-axis]: (0.000,0.000,1.000)
            Rel. specific ang. mom. (h):         1.0000

                       Eccentricity (e):         0.8599
                  Semi-latus rectum (p):         1.8599
                     Semimajor axis (a):         7.1388
                     Semiminor axis (b):         3.6438
       Orbital inclination (psi_0; deg):         0.0000
      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (0.860,0.000,0.000)
             Apsidal angle (phi_0; deg):         0.0000
   Argument of periapsis (omega_0; deg):         0.0000
                   Rel. pericentre (rp):         1.0000
                Vel. at pericentre (vp):         1.0000
                    Rel. apocentre (ra):        13.2775
                 Vel. at apocentre (va):         0.0753
                Rel. orbital period (T):       163.4408
        Approx. orbit circumference (u):        34.7642
        Approx. pericentric period (Tp):        34.7642
     Rel. specific potential energy (V):    -2.9569E+00
       Rel. specific kinetic energy (T):     5.0000E-01
         Rel. specific total energy (E):    -2.4569E+00


Orbit integration (forward):
	Time range [t0,t1] = [0.0,10.0]
	Time steps 10000;  Step size: 1.00E-03
	Maximum time step to avoid energy drift: 7.824692E+00
	Expected accumulated error in x and v of order 1.00E-04

WARNING:
	longitude of ascending node undefined. Set to 0 by convention.


Output with timestep frequency output_freq=10 written to file:

	./output/params_ini_3d_two_body_22_out.dat

Conservation laws:
                          Energy conservation to better than 1.265E-06 %.
    Angular momentum conservation (magnitude) to better than 4.445E-11 %.
    Angular momentum conservation (direction) to better than 0.000E+00 %.

A series of relevant plots for this system can be produced using:

shell> gnuplot -e 'projPlane = "op"; timeStep = 1.e-3; timeFreq = 10; ffw=1; plotRelOrbit = "T"; dataFile="./output/params_ini_3d_two_body_22_out.dat"' plot_two_body_orbit.gp

Try also setting plotRelOrbit = "F".


'''

from utils import funcs


t_0 = 0.									# initial time (Gyr)
t_1 = 10.									# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3							# integration time step

# Body 1
Mass1_scale = 1.0e6							# mass scale (Msun)
rs1 = 2										# Homogeneous sphere radius
Potential1 = funcs.Sphere_Potential(Mass1_scale,rs1)	# potential (km/s)^2
x1_0 = 0.									# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.									# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1									# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)	# potential (km/s)^2
x2_0 = 1.											# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.											# velocities (km/s):
vy2_0 = 1.
vz2_0 = 0.
