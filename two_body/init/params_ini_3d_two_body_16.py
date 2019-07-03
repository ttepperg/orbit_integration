'''
Author:	Thorsten Tepper Garcia
Date:	26/06/2019

This set of parameters produces a well-behaved two-body system, with a stable, highly
precessing orbit characterised by the following parameters:

                           Reduced mass:     9.9999E+00

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.0000
      Initial tangential vel. (v_tan_0):         1.0000
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,1.000)
                           [normalised]: (0.000,0.000,1.000)
            Rel. specific ang. mom. (h):         1.0000

             Rel. semi-latus rectum (p):         0.1916
       Rel. eccentricity vector (e_vec): (-0.808,0.000,0.000)
                  Rel. eccentricity (e):         0.8084
                Rel. semimajor axis (a):         0.5530
                Rel. semiminor axis (b):         0.3255
                   Rel. pericentre (rp):         0.1059
                Vel. at pericentre (vp):         9.4391
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         1.0000
                Rel. orbital period (T):         1.1309
        Approx. orbit circumference (u):         3.5327
        Approx. pericentric period (Tp):         0.3743
         Approx. apocentric period (Ta):         3.5327
             Apsidal angle (phi_0; deg):       180.0000
       Orbital inclination (psi_0; deg):         0.0000

              Rel. potential energy (T):    -1.8731E+02
                Rel. kinetic energy (T):     5.0000E+00
                  Rel. total energy (T):    -1.8231E+02


           Energy conservation to better than 1.343E-04 %.
 Angular momentum conservation to better than 2.284E-10 %.
 

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
Mass1 = Mass1_cum(Rvir1)										# total ('virial') mass (Msun)
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
vy2_0 = 1.
vz2_0 = 0.
