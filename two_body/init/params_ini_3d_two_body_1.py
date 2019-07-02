'''
Author:	Thorsten Tepper Garcia
Date:	22/06/2019

This set of parameters produces a well-behaved equal mass, two-body system, with a stable, non-precessing, inclined, elliptical orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+01
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  5.0000E+00
      Initial rel. distance (r210):      0.0200
         Initial rel. speed (v210):      0.0566
  Initial tangential vel. (v0_tan):      0.0566
      Initial radial vel. (v0_rad):      0.0000
              True anomaly (theta):    180.0000
  Rel. specific ang. mom. vec. (h): (-0.000,0.001,0.001)
       Rel. specific ang. mom. (h):      0.0011
              Rel. mean radius (p):      0.0149
             Rel. eccentricity (e):      0.2560
           Rel. semimajor axis (a):      0.0159
           Rel. semiminor axis (b):      0.0154
              Rel. pericentre (rp):      0.0118
               Rel. apocentre (ra):      0.0200
           Rel. orbital period (T):      1.3613
         Rel. potential energy (T): -2.1505E-02
           Rel. kinetic energy (T):  8.0000E-03
             Rel. total energy (T): -1.3505E-02


           Energy conservation to better than 8.418E-04 %.
 Angular momentum conservation to better than 2.449E-07 %.


A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "F"; projPlane = "xy"' plot_orbit_two_body.gp

Try also projPlane = "yz"

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
Mass1 = 1.0e1													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(amp=pc.Grav*Mass1)			# potential (km/s)^2
x1_0 = -0.01													# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = -0.02
vz1_0 = 0.02

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = 0.01														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 0.02
vz2_0 = -0.02
