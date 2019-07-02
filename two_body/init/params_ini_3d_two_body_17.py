'''
Author:	Thorsten Tepper Garcia
Date:	26/06/2019

This set of parameters produces a well-behaved two-body system, with a stable, nearly perfect
circular orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+06
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  9.9999E+00
     Initial rel. distance (r21_0):      1.0000
        Initial rel. speed (v21_0):      2.0739
 Initial tangential vel. (v_tan_0):      2.0739
     Initial radial vel. (v_rad_0):      0.0000
Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,2.074)
       Rel. specific ang. mom. (h):      2.0739
              Rel. mean radius (p):      1.0000
  Rel. eccentricity vector (e_vec): (-0.000,0.000,0.000)
             Rel. eccentricity (e):      0.0000
           Rel. semimajor axis (a):      1.0000
           Rel. semiminor axis (b):      1.0000
              Rel. pericentre (rp):      1.0000
               Rel. apocentre (ra):      1.0000
 Orbital rotation angle (thetha_0):    180.0000
           Rel. orbital period (T):      3.0296
         Rel. potential energy (T): -4.3009E+01
           Rel. kinetic energy (T):  2.1504E+01
             Rel. total energy (T): -2.1505E+01


           Energy conservation to better than 3.873E-09 %.
 Angular momentum conservation to better than 1.859E-11 %.

A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "T"; projPlane = "xy"; velVectorScale=0.1' plot_orbit_two_body.gp

Use the following settings in plot_orbit_two_body.gp

lengthUnitName = "kpc"
lengthUnit = 1.
velUnitName = "kms^{-1}"
velUnit = 1.
timeUnitName = "Gyr"
timeUnit = 0.978
timeFreq = 10
timeStep = 0.001
pauseStep = 0.001

'''

import config.phys_consts as pc
from utils import funcs


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e1														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Body 1
Mass1 = 1.e06													# total mass (Msun)
rs1=2.e-1
Potential1 = funcs.Hernquist_Potential(amp=pc.Grav*Mass1,a=rs1)
Mass1_cum =  funcs.Hernquist_Mass(mass=Mass1,a=rs1)				# cumulative mass
x1_0 = 0.
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2)							# 'cumulative' mass, trivially equal to Mass2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = funcs.v_circ_gen(x2_0-x1_0,y2_0-y1_0,z2_0-z1_0,pot=Potential1)
vz2_0 = 0.
