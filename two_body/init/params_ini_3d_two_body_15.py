'''
Author:	Thorsten Tepper Garcia
Date:	26/06/2019

N.B.: This is the time reversed evolution of the system defined by params_ini_3d_two_body_13.py.

This set of parameters produces a well-behaved two-body system, with a stable, highly precessing
orbit characterised by the following parameters:

              Total mass of body 1:  1.2136E+06
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  9.9999E+00
     Initial rel. distance (r21_0):      0.8032
        Initial rel. speed (v21_0):      1.7902
 Initial tangential vel. (v_tan_0):      1.2441
     Initial radial vel. (v_rad_0):      1.2874
Rel. specific ang. mom. vec. (h_vec): (-0.000,-0.000,-0.999)
       Rel. specific ang. mom. (h):      0.9992
              Rel. mean radius (p):      0.2563
  Rel. eccentricity vector (e_vec): (-0.685,0.322,0.000)
             Rel. eccentricity (e):      0.7567
           Rel. semimajor axis (a):      0.5998
           Rel. semiminor axis (b):      0.3921
              Rel. pericentre (rp):      0.1459
               Rel. apocentre (ra):      1.0537
 Orbital rotation angle (thetha_0):    154.8044
           Rel. orbital period (T):      1.4789
         Rel. potential energy (T): -1.9836E+02
           Rel. kinetic energy (T):  1.6025E+01
             Rel. total energy (T): -1.8233E+02


           Energy conservation to better than 1.371E-04 %.
 Angular momentum conservation to better than 2.152E-10 %.

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
c1 = 1.0e0														# NFW concentration
rs1 = 1.0e0														# NFW scale radius (kpc)
rho01 = 5.0e5													# core density (Msun/kpc**3)
Rvir1 = c1*rs1
Potential1 = funcs.NFW_Potential(rho01,rs1)						# potential (km/s)^2
Mass1_cum = funcs.NFW_Mass(rho01,rs1)
Mass1 = Mass1_cum(Rvir1)											# total ('virial') mass (Msun)
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2)							# 'cumulative' mass, trivially equal to Mass2
x2_0 = 5.0467E-01												# positions (kpc)
y2_0 = -6.2483E-01
z2_0 = 8.9603E-13
vx2_0 = -1.7767E+00												# velocities (km/s):
vy2_0 = 2.1981E-01
vz2_0 = 1.3800E-12
