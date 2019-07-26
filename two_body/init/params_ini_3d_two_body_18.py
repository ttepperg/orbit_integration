'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable,
non-precessing inclined, elliptical orbit characterised by the following parameters:

                           Reduced mass:     5.0000E+05

          Initial rel. distance (r21_0):         1.7321
             Initial rel. speed (v21_0):         2.0000
      Initial tangential vel. (v_tan_0):         1.6330
          Initial radial vel. (v_rad_0):         1.1547
   Rel. specific ang. mom. vec. (h_vec): (-2.000,0.000,2.000)
                           [normalised]: (-0.707,0.000,0.707)
            Rel. specific ang. mom. (h):         2.8284

             Rel. semi-latus rectum (p):         0.9300
       Rel. eccentricity vector (e_vec): (-0.112,-0.577,-0.112)
                  Rel. eccentricity (e):         0.5988
                Rel. semimajor axis (a):         1.4499
                Rel. semiminor axis (b):         1.1612
                   Rel. pericentre (rp):         0.5817
                Vel. at pericentre (vp):         4.8623
                    Rel. apocentre (ra):         2.3182
                 Vel. at apocentre (va):         1.2201
                Rel. orbital period (T):         3.7403
        Approx. orbit circumference (u):         9.1381
        Approx. pericentric period (Tp):         1.8794
         Approx. apocentric period (Ta):         7.4896
             Apsidal angle (phi_0; deg):      -105.3847
       Orbital inclination (psi_0; deg):        45.0000

              Rel. potential energy (T):    -2.4831E+06
                Rel. kinetic energy (T):     1.0000E+06
                  Rel. total energy (T):    -1.4831E+06


           Energy conservation to better than 2.003E-03 %.
 Angular momentum conservation to better than 2.367E-10 %.

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
Mass1_scale = 1.0e6													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(mass=Mass1_scale)					# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = -1.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e6													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(mass=Mass2_scale)					# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 1.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 2.
vz2_0 = 0.
