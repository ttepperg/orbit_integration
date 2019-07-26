'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, nearly perfect
circular orbit characterised by the following parameters:

                   Total mass of body 1:     6.2832E+06
                   Total mass of body 2:     1.0000E+01
                           Reduced mass:     1.0000E+01

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         2.2846
      Initial tangential vel. (v_tan_0):         2.2846
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,-2.285)
                           [normalised]: (0.000,0.000,-1.000)
                         [along z-axis]: (0.000,0.000,-2.285)
            Rel. specific ang. mom. (h):         2.2846

                       Eccentricity (e):         0.0000
                  Semi-latus rectum (p):         1.0000
                     Semimajor axis (a):         1.0000
                     Semiminor axis (b):         1.0000
       Orbital inclination (psi_0; deg):       180.0000
                 Ascending node (n_vec): (-0.000,0.000,0.000)

WARNING: longitude of ascending node undefined. Set to 0 by convention.
      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (0.000,0.000,0.000)
             Apsidal angle (phi_0; deg):         0.0000
   Argument of periapsis (omega_0; deg):         0.0000
                   Rel. pericentre (rp):         1.0000
                Vel. at pericentre (vp):         2.2846
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         2.2846
                Rel. orbital period (T):         2.7502
        Approx. orbit circumference (u):         6.2832
        Approx. pericentric period (Tp):         2.7502
     Rel. specific potential energy (V):    -1.8731E+01
       Rel. specific kinetic energy (T):     2.6098E+00
         Rel. specific total energy (E):    -1.6122E+01


Conservation laws:
                          Energy conservation to better than 5.254E-11 %.
    Angular momentum conservation (magnitude) to better than 8.442E-11 %.
    Angular momentum conservation (direction) to better than 0.000E+00 %.


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

from utils import funcs
from math import pi


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e1														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Body 1
rs1 = 1.0e0														# NFW scale radius (kpc)
rho01 = 5.0e5													# core density (Msun/kpc**3)
Mass1_scale = 4. * pi * rho01 * rs1**3
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)						# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)						# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 0.
vz2_0 = 0.


# redefine velocities to obtain a circular orbit
# assumption: centre of mass at rest
from math import sqrt
Mass1_cum = funcs.NFW_Mass(Mass1_scale,rs1)
Mass2_cum = funcs.Kepler_Mass(Mass2_scale)							# 'cumulative' mass, trivially equal to Mass2_scale
r0 = [x2_0-x1_0,y2_0-y1_0,z2_0-z1_0]
M1 = Mass1_cum(*r0)
M2 = Mass2_cum(*r0)
mu = M1*M2/(M1+M2)
gx,gy,gz = funcs.grav_field(*r0,pot=Potential1)
g = [gx,gy,gz]
vel_rel = sqrt(funcs.norm(*r0)*M2*funcs.norm(*g) / mu)
vy1_0 = vel_rel * mu / M1
vy2_0 = -1.* vel_rel * mu / M2
