'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved equal mass, two-body system, with a stable, non-precessing nearly perfect circular orbit characterised by the following parameters:

                           Reduced mass:     5.0000E+00

          Initial rel. distance (r21_0):         0.0200
             Initial rel. speed (v21_0):         0.0656
      Initial tangential vel. (v_tan_0):         0.0656
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,-0.001)
                           [normalised]: (0.000,0.000,-1.000)
            Rel. specific ang. mom. (h):         0.0013

             Rel. semi-latus rectum (p):         0.0200
       Rel. eccentricity vector (e_vec): (-0.000,0.000,0.000)
                  Rel. eccentricity (e):         0.0000
                Rel. semimajor axis (a):         0.0200
                Rel. semiminor axis (b):         0.0200
                   Rel. pericentre (rp):         0.0200
                Vel. at pericentre (vp):         0.0656
                    Rel. apocentre (ra):         0.0200
                 Vel. at apocentre (va):         0.0656
                Rel. orbital period (T):         1.9161
        Approx. orbit circumference (u):         0.1257
        Approx. pericentric period (Tp):         1.9161
         Approx. apocentric period (Ta):         1.9161
             Apsidal angle (phi_0; deg):       180.0000
       Orbital inclination (psi_0; deg):       180.0000

              Rel. potential energy (T):    -2.1505E-02
                Rel. kinetic energy (T):     1.0752E-02
                  Rel. total energy (T):    -1.0752E-02


           Energy conservation to better than 5.715E-08 %.
 Angular momentum conservation to better than 2.734E-08 %.


A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "F"; projPlane = "xy"; velVectorScale=0.1' plot_orbit_two_body.gp

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
delta_t = 1.0e-3												# time step

# Body 1
Mass1_scale = 1.0e1													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(mass=Mass1_scale)					# potential (km/s)^2
x1_0 = -0.01													# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# arbitrary velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2_scale)					# potential (km/s)^2
x2_0 = 0.01														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# arbitrary velocities (km/s):
vy2_0 = 0.
vz2_0 = 0.


# redefine velocities to obtain a circular orbit
# assumption: centre of mass at rest
from math import sqrt
M1 = Mass1_scale
M2 = Mass2_scale
mu = M1*M2/(M1+M2)
r0 = [x2_0-x1_0,y2_0-y1_0,z2_0-z1_0]
gx,gy,gz = funcs.grav_field(*r0,pot=Potential1)
g = [gx,gy,gz]
vel_rel = sqrt(funcs.norm(*r0)*M2*funcs.norm(*g) / mu)
vy1_0 = vel_rel * mu / M1
vy2_0 = -1.* vel_rel * mu / M2
