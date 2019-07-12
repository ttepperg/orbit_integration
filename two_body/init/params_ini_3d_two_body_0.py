'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved equal mass, two-body system, with a stable,
non-precessing elliptical orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+01
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  5.0000E+00
     Initial rel. distance (r21_0):      0.0200
        Initial rel. speed (v21_0):      0.0400
 Initial tangential vel. (v_tan_0):      0.0400
     Initial radial vel. (v_rad_0):      0.0000
Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,0.001)
       Rel. specific ang. mom. (h):      0.0008
              Rel. mean radius (p):      0.0074
  Rel. eccentricity vector (e_vec): (-0.628,0.000,0.000)
             Rel. eccentricity (e):      0.6280
           Rel. semimajor axis (a):      0.0123
           Rel. semiminor axis (b):      0.0096
              Rel. pericentre (rp):      0.0046
               Rel. apocentre (ra):      0.0200
 Orbital rotation angle (thetha_0):    180.0000
           Rel. orbital period (T):      0.9225
         Rel. potential energy (T): -2.1505E-02
           Rel. kinetic energy (T):  4.0000E-03
             Rel. total energy (T): -1.7505E-02

           Energy conservation to better than 4.548E-02 %.
 Angular momentum conservation to better than 3.300E-06 %.


A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "F"; projPlane = "xy"' plot_orbit_two_body.gp

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
Potential1 = funcs.Kepler_Potential(mass=Mass1)					# potential (km/s)^2
Mass1_cum = funcs.Kepler_Mass(Mass1)							# 'cumulative' mass, trivially equal to Mass1
x1_0 = -0.01													# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = -0.02
vz1_0 = 0.

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2)							# 'cumulative' mass, trivially equal to Mass1
x2_0 = 0.01														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 0.02
vz2_0 = 0.
