'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, non-precessing
orbit characterised by the following parameters:

                    Potential of body 1:     Kepler_Pot
                    Potential of body 2:     Kepler_Pot
                   Mass scale of body 1:     1.0000E+06
                   Mass scale of body 2:     1.0000E+01

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.5677
      Initial tangential vel. (v_tan_0):         1.5677
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,1.568)
                           [normalised]: (0.000,0.000,1.000)
                         [along z-axis]: (0.000,0.000,1.568)
            Rel. specific ang. mom. (h):         1.5677

                       Eccentricity (e):         0.4286
                  Semi-latus rectum (p):         0.5714
                     Semimajor axis (a):         0.7000
                     Semiminor axis (b):         0.6325
       Orbital inclination (psi_0; deg):         0.0000
                 Ascending node (n_vec): (0.000,0.000,0.000)

WARNING:
	longitude of ascending node undefined. Set to 0 by convention.

      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (-0.429,0.000,0.000)
             Apsidal angle (phi_0; deg):       180.0000
   Argument of periapsis (omega_0; deg):       180.0000
                   Rel. pericentre (rp):         0.4000
                Vel. at pericentre (vp):         3.9193
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         1.5677
                Rel. orbital period (T):         1.7744
        Approx. orbit circumference (u):         4.1887
        Approx. pericentric period (Tp):         1.0688
     Rel. specific potential energy (V):    -4.3010E+00
       Rel. specific kinetic energy (T):     1.2288E+00
         Rel. specific total energy (E):    -3.0721E+00


Conservation laws:
                          Energy conservation to better than 1.816E-03 %.
    Angular momentum conservation (magnitude) to better than 7.685E-11 %.
    Angular momentum conservation (direction) to better than 0.000E+00 %.


A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "T"; projPlane = "xy"' plot_orbit_two_body.gp

Use the following settings in plot_orbit_two_body.gp

lengthUnitName = "kpc"
lengthUnit = 1.
velUnitName = "kms^{-1}"
velUnit = 1.
velVectorScale = 0.1	# for nice arrows on plot
timeUnitName = "Gyr"
timeUnit = 0.978
timeFreq = 10
timeStep = 0.001
pauseStep = 0.001

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.														# initial time (Gyr)
t_1 = 10.														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Body 1
Mass1_scale = 1.0e6													# mass scale (Msun)
Potential1 = funcs.Kepler_Potential(Mass1_scale)					# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)					# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
sm_a = 0.7														# must be > r_0/2
vy2_0 = \
	funcs.vis_viva(x2_0,y2_0,z2_0,amp=pc.Grav*(Mass1_scale+Mass2_scale),a=sm_a)	# elliptical orbit with semi-major axis a > r/2
vz2_0 = 0.
