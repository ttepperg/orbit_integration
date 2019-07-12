'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable,
hyperbolic (nearly parabolic though) orbit characterised by the following parameters:

              Total mass of body 1:  1.0000E+06
              Total mass of body 2:  1.0000E+01
                      Reduced mass:  9.9999E+00
     Initial rel. distance (r21_0):      2.0000
        Initial rel. speed (v21_0):      2.2361
 Initial tangential vel. (v_tan_0):      1.0000
     Initial radial vel. (v_rad_0):      2.0000
Rel. specific ang. mom. vec. (h_vec): (-0.000,0.000,2.000)
       Rel. specific ang. mom. (h):      2.0000
              Rel. mean radius (p):      0.9300
  Rel. eccentricity vector (e_vec): (0.930,0.535,0.000)
             Rel. eccentricity (e):      1.0729
           Rel. semimajor axis (a):     -6.1527
           Rel. semiminor axis (b):      2.3921
              Rel. pericentre (rp):      0.4487
               Rel. apocentre (ra):         nan
 Orbital rotation angle (thetha_0):     29.9092
           Rel. orbital period (T):         nan
         Rel. potential energy (T): -2.1505E+01
           Rel. kinetic energy (T):  2.5000E+01
             Rel. total energy (T):  3.4952E+00


           Energy conservation to better than 1.722E-02 %.
 Angular momentum conservation to better than 7.544E-11 %.

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
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e1													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2)					# potential (km/s)^2
x2_0 = 0.														# positions (kpc)
y2_0 = -2.
z2_0 = 0.
vx2_0 = 1.														# velocities (km/s):
vy2_0 = 2.
vz2_0 = 0.
