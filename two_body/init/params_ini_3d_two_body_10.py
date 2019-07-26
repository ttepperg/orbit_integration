'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, non-precessing
orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+06
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  9.9999E+00
     Initial rel. distance (r21_0):      1.0050
        Initial rel. speed (v21_0):      2.2361
 Initial tangential vel. (v_tan_0):      1.1940
     Initial radial vel. (v_rad_0):      1.8906
  Rel. specific ang. mom. vec. (h): (-0.000,0.000,1.200)
       Rel. specific ang. mom. (h):      1.2000
              Rel. mean radius (p):      0.3348
  Rel. eccentricity vector (e_vec): (0.459,0.716,0.000)
             Rel. eccentricity (e):      0.8503
           Rel. semimajor axis (a):      1.2084
           Rel. semiminor axis (b):      0.6361
              Rel. pericentre (rp):      0.1810
               Rel. apocentre (ra):      2.2358
 Orbital rotation angle (thetha_0):     57.3664
           Rel. orbital period (T):      4.0245
         Rel. potential energy (T): -4.2796E+01
           Rel. kinetic energy (T):  2.5000E+01
             Rel. total energy (T): -1.7796E+01


           Energy conservation to better than 1.090E-01 %.
 Angular momentum conservation to better than 4.026E-11 %.

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

from utils import funcs


t_0 = 0.														# initial time (Gyr)
t_1 = 10.														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Body 1
Mass1_scale = 1.0e6													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(mass=Mass1_scale)					# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2_scale)					# potential (km/s)^2
x2_0 = 0.1														# positions (kpc)
y2_0 = -1.
z2_0 = 0.
vx2_0 = 1.														# velocities (km/s):
vy2_0 = 2.
vz2_0 = 0.
