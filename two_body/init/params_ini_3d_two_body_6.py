'''
Author:	Thorsten Tepper Garcia
Date:	22/06/2019

This set of parameters produces a well-behaved equal mass, two-body system, with a stable, inclined,
higly precessing orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+06
              Total mass of body 2:  1.0000E+06
                      Reduced mass:  5.0000E+05
      Initial rel. distance (r210):      1.4142
         Initial rel. speed (v210):      1.0000
  Initial tangential vel. (v0_tan):      1.0000
      Initial radial vel. (v0_rad):      0.0000
              True anomaly (theta):    180.0000
  Rel. specific ang. mom. vec. (h): (-1.000,0.000,1.000)
       Rel. specific ang. mom. (h):      1.4142
              Rel. mean radius (p):      0.2325
             Rel. eccentricity (e):      0.8356
           Rel. semimajor axis (a):      0.7704
           Rel. semiminor axis (b):      0.4232
              Rel. pericentre (rp):      0.1267
               Rel. apocentre (ra):      1.4142
           Rel. orbital period (T):      1.4487
         Rel. potential energy (T): -3.0374E+06
           Rel. kinetic energy (T):  2.5000E+05
             Rel. total energy (T): -2.7874E+06

Time range [t0,t1] = [0.0,10.0]
Time steps 10000;  Step size: 0.001

Writing output to file ./output/two_body_orbit_int_gen_3d_leap_phys.dat with timestep frequency 10


           Energy conservation to better than 2.591E-01 %.
 Angular momentum conservation to better than 5.007E-11 %.


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
rs = 1.e-1														# Plummer scale radius
Potential1 = funcs.Plummer_Potential(amp=pc.Grav*Mass1,a=rs)	# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e6													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 1.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 1.
vz2_0 = 0.
