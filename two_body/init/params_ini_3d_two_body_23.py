'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved equal mass, two-body system, with a point mass and an extended body on a precessing, elliptic orbit, characterised by the following parameters:


          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.0000
      Initial tangential vel. (v_tan_0):         1.0000
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,1.000)
                           [normalised]: (0.000,0.000,1.000)
                         [along z-axis]: (0.000,0.000,1.000)
            Rel. specific ang. mom. (h):         1.0000

                       Eccentricity (e):         0.6751
                  Semi-latus rectum (p):         0.3249
                     Semimajor axis (a):         0.5970
                     Semiminor axis (b):         0.4404
       Orbital inclination (psi_0; deg):         0.0000
                 Ascending node (n_vec): (0.000,0.000,0.000)
      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (-0.675,0.000,0.000)
             Apsidal angle (phi_0; deg):       180.0000
   Argument of periapsis (omega_0; deg):       180.0000
                   Rel. pericentre (rp):         0.1940
                Vel. at pericentre (vp):         5.1551
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         1.0000
                Rel. orbital period (T):         1.6521
        Approx. orbit circumference (u):         3.2777
        Approx. pericentric period (Tp):         0.6358
     Rel. specific potential energy (V):    -3.8469E+00
       Rel. specific kinetic energy (T):     5.0000E-01
         Rel. specific total energy (E):    -3.3469E+00

Conservation laws:
                          Energy conservation to better than 2.574E-04 %.
    Angular momentum conservation (magnitude) to better than 5.751E-11 %.
    Angular momentum conservation (direction) to better than 0.000E+00 %.

A series of relevant plots for this system can be produced using:

shell> gnuplot -e 'projPlane = "op"; timeStep = 1.e-3; timeFreq = 10; ffw=1; plotRelOrbit = "T"; dataFile="./output/params_ini_3d_two_body_23_out.dat"' plot_two_body_orbit.gp

Try also setting plotRelOrbit = "F".


'''

from utils import funcs


t_0 = 0.									# initial time (Gyr)
t_1 = 10.									# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3							# integration time step

# Body 1
Mass1_scale = 1.0e6							# mass scale (Msun)
rs1 = 0.5									# Plummer scale radius
Potential1 = funcs.Plummer_Potential(Mass1_scale,rs1)	# potential (km/s)^2
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
