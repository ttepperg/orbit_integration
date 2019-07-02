'''
Author:	Thorsten Tepper Garcia
Date:	26/06/2019

This set of parameters produces a well-behaved two-body system, with a stable, non-precessing
elliptical orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+06
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  9.9999E+00
     Initial rel. distance (r21_0):      1.0000
        Initial rel. speed (v21_0):      1.0000
 Initial tangential vel. (v_tan_0):      1.0000
     Initial radial vel. (v_rad_0):      0.0000
Rel. specific ang. mom. vec. (h_vec): (0.000,-1.000,0.000)
       Rel. specific ang. mom. (h):      1.0000
              Rel. mean radius (p):      0.2325
  Rel. eccentricity vector (e_vec): (-0.767,0.000,-0.000)
             Rel. eccentricity (e):      0.7675
           Rel. semimajor axis (a):      0.5658
           Rel. semiminor axis (b):      0.3627
              Rel. pericentre (rp):      0.1315
               Rel. apocentre (ra):      1.0000
 Orbital rotation angle (thetha_0):    180.0000
           Rel. orbital period (T):      1.2893
         Rel. potential energy (T): -4.3009E+01
           Rel. kinetic energy (T):  5.0000E+00
             Rel. total energy (T): -3.8009E+01


           Energy conservation to better than 1.715E-01 %.
 Angular momentum conservation to better than 1.076E-10 %.


A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "T"; projPlane = "xz"' plot_orbit_two_body.gp

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
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 0.
vz2_0 = 1.
