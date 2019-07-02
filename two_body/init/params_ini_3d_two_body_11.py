'''
Author:	Thorsten Tepper Garcia
Date:	22/062019

This set of parameters produces a well-behaved two-body system, with a stable,
hyperbolic orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+06
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  9.9999E+00
     Initial rel. distance (r21_0):      2.0000
        Initial rel. speed (v21_0):      2.8284
 Initial tangential vel. (v_tan_0):      2.0000
     Initial radial vel. (v_rad_0):      2.0000
Rel. specific ang. mom. vec. (h_vec): (-0.000,0.000,4.000)
       Rel. specific ang. mom. (h):      4.0000
              Rel. mean radius (p):      3.7201
  Rel. eccentricity vector (e_vec): (1.860,-0.860,0.000)
             Rel. eccentricity (e):      2.0493
           Rel. semimajor axis (a):     -1.1627
           Rel. semiminor axis (b):      2.0798
              Rel. pericentre (rp):      1.2200
               Rel. apocentre (ra):         nan
 Orbital rotation angle (thetha_0):    -24.8149
           Rel. orbital period (T):         nan
         Rel. potential energy (T): -2.1505E+01
           Rel. kinetic energy (T):  4.0000E+01
             Rel. total energy (T):  1.8495E+01


           Energy conservation to better than 1.082E-04 %.
 Angular momentum conservation to better than 2.478E-11 %.

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
Mass1 = 1.0e6													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(amp=pc.Grav*Mass1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = 0.														# positions (kpc)
y2_0 = -2.
z2_0 = 0.
vx2_0 = 2.														# velocities (km/s):
vy2_0 = 2.
vz2_0 = 0.
