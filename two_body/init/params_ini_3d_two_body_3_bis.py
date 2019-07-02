'''
Author:	Thorsten Tepper Garcia
Date:	22/06/2019

This set of parameters produces a well-behaved two-body system, with a stable, non-precessing
elliptical orbit characterised by the following parameters:

                           Reduced mass:     5.0000E+05

          Initial rel. distance (r21_0):         1.4142
             Initial rel. speed (v21_0):         2.0000
      Initial tangential vel. (v_tan_0):         1.4142
          Initial radial vel. (v_rad_0):         1.4142
   Rel. specific ang. mom. vec. (h_vec): (0.000,0.000,2.000)
                           [normalised]: (0.000,0.000,1.000)
            Rel. specific ang. mom. (h):         2.0000

             Rel. semi-latus rectum (p):         0.4650
       Rel. eccentricity vector (e_vec): (-0.242,-0.707,0.000)
                  Rel. eccentricity (e):         0.7474
                Rel. semimajor axis (a):         1.0535
                Rel. semiminor axis (b):         0.6999
                   Rel. pericentre (rp):         0.2661
                Vel. at pericentre (vp):         7.5154
                    Rel. apocentre (ra):         1.8409
                 Vel. at apocentre (va):         1.0864
                Rel. orbital period (T):         2.3166
        Approx. orbit circumference (u):         6.6868
        Approx. pericentric period (Tp):         0.8897
         Approx. apocentric period (Ta):         6.1549
             Apsidal angle (phi_0; deg):      -108.8995
       Orbital inclination (psi_0; deg):         0.0000

              Rel. potential energy (T):    -3.0412E+06
                Rel. kinetic energy (T):     1.0000E+06
                  Rel. total energy (T):    -2.0412E+06


           Energy conservation to better than 3.759E-02 %.
 Angular momentum conservation to better than 3.361E-11 %.


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

import config.phys_consts as pc
from utils import funcs


t_0 = 0.														# initial time (Gyr)
t_1 = 10.														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-3												# integration time step

# Body 1
Mass1 = 1.0e6													# total mass (Msun)
Potential1 = funcs.Kepler_Potential(amp=pc.Grav*Mass1)			# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.

# Body 2
Mass2 = 1.0e6													# total mass (Msun)
Potential2 = funcs.Kepler_Potential(amp=pc.Grav*Mass2)			# potential (km/s)^2
x2_0 = 1.														# positions (kpc)
y2_0 = 1.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 2.
vz2_0 = 0.
