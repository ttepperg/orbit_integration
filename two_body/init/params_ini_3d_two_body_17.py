'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, nearly perfect
circular orbit characterised by the following parameters:

                           Reduced mass:     9.9999E+00

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.7282
      Initial tangential vel. (v_tan_0):         1.7282
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,-1.728)
                           [normalised]: (0.000,0.000,-1.000)
                         [along z-axis]: (0.000,0.000,-1.728)
            Rel. specific ang. mom. (h):         1.7282

                       Eccentricity (e):         0.0000
                  Semi-latus rectum (p):         1.0000
                     Semimajor axis (a):         1.0000
                     Semiminor axis (b):         1.0000
       Orbital inclination (psi_0; deg):       180.0000
                 Ascending node (n_vec): (-0.000,0.000,0.000)
      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (-0.000,0.000,0.000)
             Apsidal angle (phi_0; deg):         0.0000
   Argument of periapsis (omega_0; deg):         0.0000
                   Rel. pericentre (rp):         1.0000
                Vel. at pericentre (vp):         1.7282
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         1.7282
                Rel. orbital period (T):         3.6356
        Approx. orbit circumference (u):         6.2832
        Approx. pericentric period (Tp):         3.6356
              Rel. potential energy (V):    -3.5841E+01
                Rel. kinetic energy (T):     1.4934E+01
                  Rel. total energy (E):    -2.0907E+01


           Energy conservation to better than 1.492E-10 %.
 Angular momentum conservation to better than 3.548E-11 %.

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


t_0 = 0.0e0														# initial time (Gyr)
t_1 = 1.0e1														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Body 1
Mass1_scale = 1.e06													# mass scale (Msun)
rs1=2.e-1
Potential1 = funcs.Hernquist_Potential(Mass1_scale,rs1)		# potential (km/s)^2
Mass1_cum =  funcs.Hernquist_Mass(Mass1_scale,rs1)				# cumulative mass
x1_0 = 0.
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)			# potential (km/s)^2
Mass2_cum = funcs.Kepler_Mass(Mass2_scale)							# 'cumulative' mass, trivially equal to Mass2_scale
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 0.
vz2_0 = 0.


# redefine velocities to obtain a circular orbit
# assumption: centre of mass at rest
from math import sqrt
r0 = [x2_0-x1_0,y2_0-y1_0,z2_0-z1_0]
M1 = Mass1_cum(*r0)
M2 = Mass2_cum(*r0)
mu = M1*M2/(M1+M2)
gx,gy,gz = funcs.grav_field(*r0,pot=Potential1)
g = [gx,gy,gz]
vel_rel = sqrt(funcs.norm(*r0)*M2*funcs.norm(*g) / mu)
vy1_0 = vel_rel * mu / M1
vy2_0 = -1.* vel_rel * mu / M2
