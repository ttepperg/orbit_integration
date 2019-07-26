'''
Author:	Thorsten Tepper Garcia

This set of parameters produces a well-behaved two-body system, with a stable, non-precessing
orbit characterised by the following parameters:

                           Reduced mass:     9.9999E+00

          Initial rel. distance (r21_0):         1.0000
             Initial rel. speed (v21_0):         1.0000
      Initial tangential vel. (v_tan_0):         1.0000
          Initial radial vel. (v_rad_0):         0.0000
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,1.000)
                           [normalised]: (0.000,0.000,1.000)
                         [along z-axis]: (0.000,0.000,1.000)
            Rel. specific ang. mom. (h):         1.0000

                       Eccentricity (e):         0.7675
                  Semi-latus rectum (p):         0.2325
                     Semimajor axis (a):         0.5658
                     Semiminor axis (b):         0.3627
       Orbital inclination (psi_0; deg):         0.0000
                 Ascending node (n_vec): (0.000,0.000,0.000)
      Long. of asc. node (Omega_0; deg):         0.0000
            Eccentricity vector (e_vec): (-0.767,0.000,0.000)
             Apsidal angle (phi_0; deg):       180.0000
   Argument of periapsis (omega_0; deg):       180.0000
                   Rel. pericentre (rp):         0.1315
                Vel. at pericentre (vp):         7.6019
                    Rel. apocentre (ra):         1.0000
                 Vel. at apocentre (va):         1.0000
                Rel. orbital period (T):         1.2893
        Approx. orbit circumference (u):         2.9517
        Approx. pericentric period (Tp):         0.3883
              Rel. potential energy (V):    -4.3009E+01
                Rel. kinetic energy (T):     5.0000E+00
                  Rel. total energy (E):    -3.8009E+01


           Energy conservation to better than 1.715E-01 %.
 Angular momentum conservation to better than 1.076E-10 %.


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
Mass1_scale = 1.0e6													# mass scale (Msun)
Potential1 = funcs.Kepler_Potential(Mass1_scale)					# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2_scale = 1.0e1													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)					# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 1.
vz2_0 = 0.
