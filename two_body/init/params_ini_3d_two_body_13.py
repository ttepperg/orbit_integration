'''
Author:	Thorsten Tepper Garcia

N.B.: This is the time reversed evolution of the system defined by params_ini_3d_two_body_15.py.

This set of parameters produces a well-behaved two-body system, with a stable, highly precessing
orbit characterised by the following parameters:

                   Total mass of body 1:     6.2832E+06
                   Total mass of body 2:     1.0000E+01
                           Reduced mass:     1.0000E+01

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.0000
      Initial tangential vel. (v_tan_0):         1.0000
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,1.000)
                           [normalised]: (0.000,0.000,1.000)
                         [along z-axis]: (0.000,0.000,1.000)
            Rel. specific ang. mom. (h):         1.0000

                       Eccentricity (e):         0.8084
                  Semi-latus rectum (p):         0.1916
                     Semimajor axis (a):         0.5530
                     Semiminor axis (b):         0.3255
       Orbital inclination (psi_0; deg):         0.0000
                 Ascending node (n_vec): (0.000,0.000,0.000)

WARNING: longitude of ascending node undefined. Set to 0 by convention.
      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (-0.808,0.000,0.000)
             Apsidal angle (phi_0; deg):       180.0000
   Argument of periapsis (omega_0; deg):       180.0000
                   Rel. pericentre (rp):         0.1059
                Vel. at pericentre (vp):         9.4391
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         1.0000
                Rel. orbital period (T):         1.1309
        Approx. orbit circumference (u):         2.8060
        Approx. pericentric period (Tp):         0.2973
     Rel. specific potential energy (V):    -1.8731E+01
       Rel. specific kinetic energy (T):     5.0000E-01
         Rel. specific total energy (E):    -1.8231E+01


Conservation laws:
                          Energy conservation to better than 1.343E-04 %.
    Angular momentum conservation (magnitude) to better than 1.780E-10 %.
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
Potential1 = funcs.NFW_Potential(Mass1_scale,rs1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2_scale)					# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 1.
vz2_0 = 0.
