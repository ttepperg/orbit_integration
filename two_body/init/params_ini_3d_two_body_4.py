'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, non-precessing
nearly circular orbit characterised by the following parameters:

                           Reduced mass:     9.9999E+00

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         2.0739
      Initial tangential vel. (v_tan_0):         2.0739
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,-2.074)
                           [normalised]: (0.000,0.000,-1.000)
            Rel. specific ang. mom. (h):         2.0739

             Rel. semi-latus rectum (p):         1.0000
       Rel. eccentricity vector (e_vec): (-0.000,0.000,0.000)
                  Rel. eccentricity (e):         0.0000
                Rel. semimajor axis (a):         1.0000
                Rel. semiminor axis (b):         1.0000
                   Rel. pericentre (rp):         1.0000
                Vel. at pericentre (vp):         2.0739
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         2.0739
                Rel. orbital period (T):         3.0297
        Approx. orbit circumference (u):         6.2832
        Approx. pericentric period (Tp):         3.0297
         Approx. apocentric period (Ta):         3.0297
             Apsidal angle (phi_0; deg):       180.0000
       Orbital inclination (psi_0; deg):       180.0000

              Rel. potential energy (T):    -4.3009E+01
                Rel. kinetic energy (T):     2.1505E+01
                  Rel. total energy (T):    -2.1505E+01


           Energy conservation to better than 4.355E-10 %.
 Angular momentum conservation to better than 4.902E-11 %.


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
Mass1 = 1.0e6													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(mass=Mass1)					# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# arbitrary velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# arbitrary velocities (km/s):
vy2_0 = 0.
vz2_0 = 0.


# redefine velocities to obtain a circular orbit
# assumption: centre of mass at rest
from math import sqrt
M1 = Mass1
M2 = Mass2
mu = M1*M2/(M1+M2)
r0 = [x2_0-x1_0,y2_0-y1_0,z2_0-z1_0]
gx,gy,gz = funcs.grav_field(*r0,pot=Potential1)
g = [gx,gy,gz]
vel_rel = sqrt(funcs.norm(*r0)*M2*funcs.norm(*g) / mu)
vy1_0 = vel_rel * mu / M1
vy2_0 = -1.* vel_rel * mu / M2
