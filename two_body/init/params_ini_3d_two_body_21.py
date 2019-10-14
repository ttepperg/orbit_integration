'''
Author:	Thorsten Tepper Garcia

This set of parameters is used to test numerical orbit precision.


A series of relevant plots for this system can be produced using:

$> gnuplot -e 'dataFile="./output/params_ini_3d_two_body_21_out.dat"; plotRelOrbit = "T"; projPlane = "op"' plot_orbit_two_body.gp

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
t_1 = 3.														# total time (time unit ~ 0.978 Gyr)
delta_t = 1.0e-4												# integration time step

# Body 1
Mass1_scale = 1.0e5													# mass scale (Msun)
Potential1 = funcs.Kepler_Potential(Mass1_scale)					# potential (km/s)^2
x1_0 = 0.														# positions (kpc)
y1_0 = 0.
z1_0 = 0.
vx1_0 = 0.														# velocities (km/s):
vy1_0 = 0.
vz1_0 = 0.4

# Body 2
Mass2_scale = 1.0e0													# mass scale (Msun)
Potential2 = funcs.Kepler_Potential(Mass2_scale)					# potential (km/s)^2
x2_0 = 0.1														# positions (kpc)
y2_0 = 0.
z2_0 = 0.
vx2_0 = 0.														# velocities (km/s):
vy2_0 = 1.
vz2_0 = 0.4
