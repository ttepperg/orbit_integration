'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, nearly
circular, slowly precessing orbit characterised by the following parameters:

                           Reduced mass:     1.0000E+01

          Initial rel. distance (r21_0):         9.2000
             Initial rel. speed (v21_0):       178.0703
      Initial tangential vel. (v_tan_0):       178.0703
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (-0.000,1638.247,0.000)
                           [normalised]: (-0.000,1.000,0.000)
                         [along z-axis]: (0.000,0.000,1638.247)
            Rel. specific ang. mom. (h):      1638.2472

                       Eccentricity (e):         0.0000
                  Semi-latus rectum (p):         9.2000
                     Semimajor axis (a):         9.2000
                     Semiminor axis (b):         9.2000
       Orbital inclination (psi_0; deg):        90.0000
                 Ascending node (n_vec): (-1638.247,-0.000,0.000)
      Long. of asc. node (Omega_0; deg):       180.0000
            Eccentricity vector (e_vec): (0.000,0.000,0.000)
             Apsidal angle (phi_0; deg):         0.0000
   Argument of periapsis (omega_0; deg):       180.0000
                   Rel. pericentre (rp):         9.2000
                Vel. at pericentre (vp):       178.0703
                    Rel. apocentre (ra):         9.2000
                 Vel. at apocentre (va):       178.0703
                Rel. orbital period (T):         0.3246
        Approx. orbit circumference (u):        57.8053
        Approx. pericentric period (Tp):         0.3246
              Rel. potential energy (V):    -4.1075E+05
                Rel. kinetic energy (T):     1.5855E+05
                  Rel. total energy (E):    -2.5220E+05


           Energy conservation to better than 3.973E-10 %.
 Angular momentum conservation to better than 3.397E-10 %.

A series of relevant plots for this system can be produced using:

$> gp -e 'dataFile="./output/two_body_orbit_int_gen_3d_leap_phys.dat"; plotRelOrbit = "T"; projPlane = "xy"; velVectorScale=0.01' plot_orbit_two_body.gp

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


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e0														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-4												# integration time step

# Body 1
Mass1_scale = 1.0e11													# total mass (Msun)
rs1 = 5.e0														# Plummer scale radius
Potential1 = funcs.Plummer_Potential(mass=Mass1_scale,a=rs1)			# potential (km/s)^2
x1_0 = -1														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2_scale)					# potential (km/s)^2
x2_0 = 8.2														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 0.
vz2_0 = 0.


# redefine velocities to obtain a circular orbit
# assumption: centre of mass at rest
from math import sqrt
Mass1_cum = funcs.Plummer_Mass(mass=Mass1_scale,a=rs1)				# cumulative mass
Mass2_cum = funcs.Kepler_Mass(Mass2_scale)							# 'cumulative' mass, trivially equal to Mass2_scale
r0 = [x2_0-x1_0,y2_0-y1_0,z2_0-z1_0]
M1 = Mass1_cum(*r0)
M2 = Mass2_cum(*r0)
mu = M1*M2/(M1+M2)
gx,gy,gz = funcs.grav_field(*r0,pot=Potential1)
g = [gx,gy,gz]
vel_rel = sqrt(funcs.norm(*r0)*M2*funcs.norm(*g) / mu)
vz1_0 = vel_rel * mu / M1
vz2_0 = -1.* vel_rel * mu / M2
