#!/opt/local/bin//python3.5

'''

Author:	Thorsten Tepper Garcia
Date:	01/07/2019


BACKGROUND

The 3D equations of motion of two bodies of mass M1 and M2 orbiting around each
others potential are:

dr1dt = v1
dv1dt = - Grad Phi2

dr2dt = v2
dv2dt = - Grad Phi1

Here,

r := position vector of body 1 or 2
v := velocity vector of body 1 or 2
Grad := gradient operator (nabla)
Phi(r) = := grav. potential of body 1 or 2

Here, we consider the case where body is possibly an extended object of 'total' mass M1,
represented by an arbitrary, static potential Phi1, and a point-like body 2 of mass M2
described by a Kepler potential (i.e. Phi1 = - G M2 / r). The latter assumption may be
relaxed in future versions of the code. Note that the 'total' mass M1 may be subject to
definition, e.g. in the case of a NFW potential for which the mass diverges with radius.

Both Phi1 and Phi2 are assumed to be central potentials that depend on the distance
bwetween the bodies only.

Even though the equation of motion is of second order (m=2) in dt,
due to its 3D vector nature, it is equivalent to a system of six
first order ODEs, for each body. Thus, the rank of the system is determined
by the order m and the dimensionality of the system D=3 as 2 * m * D = 12 in
this case.

Thus, the above vectorial equations can be written down as a system of 12
coupled ordinary differential equations of first order in dt:

dx1dt = vx1
dy1dt = vy1
dz1dt = vz1
dvx1dt = - Grad_x Phi2(|r12|)
dvy1dt = - Grad_y Phi2(|r12|)
dvz1dt = - Grad_z Phi2(|r12|)

dx2dt = vx2
dy2dt = vy2
dz2dt = vz2
dvx2dt = - Grad_x Phi1(|r12|)
dvy2dt = - Grad_y Phi1(|r12|)
dvz2dt = - Grad_z Phi1(|r12|)

Here, r12 = r1 - r2. Note that the force is radial and directed towards
the centre of the potential at all times, by assumption.

Note that a 2D system can be simulated equally well with this code by simply
setting all z-components to 0 initially.


TIME INTEGRATION

The system of equations is numerically integrated using a 2nd order
leapfrog (kick-drift-kick) method. This scheme is symplectic and thus features
desired properties such  as time reversibility, conservation of the relative total
energy to any desired accuracy -- even over significant integration times --, if the
time step is chosen small enough, and conservation of the specific angular momentum
generally to machine precision.

Beware that the scheme may lead to orbit precession if the integration step is
not small enough. As a general rule of thumb, the time step delta_t should
satisfy:

	delta_t < T / sqrt(2) pi ~ 0.225 T,

where T is the (shortest) period of the system (see
https://en.wikipedia.org/wiki/Energy_drift).

Also, orbit precession may naturally result from central potentials other than
point-like masses (Bertrand's Theorem; https://en.wikipedia.org/wiki/Bertrand%27s_theorem)
and it thus do not necessarily indicate an error in the integration scheme.


SPATIAL DERIVATIVES

A number of first order, partial, spatial derivatives are required when calculating
the force corresponding to a given potential. Two methods are implemented here:

1) A central finite difference scheme of orders 2 and 4;

2) A forward finite difference scheme of orders 2 and 3.

Method 1 is preferred (and is the default), in particular because it leads to a
better conservation of angular momentum for extrem potentials such as a
point-like mass potential. Note that using a different scheme requires
to change the invoked scheme explicitly in the calculation of the accelerations.

An example on how to invoke the differentiation routines is as follows:

1) To calculate the first partial derivative of f with respect to x (var=0) at
point r = [x,y,z], using a CFD with an integration step 1.e-4 at a precision 3**2 use:
fwd_diff_first(*r, var = 0, func = f, delta_x = 1.e-4, order = 3)

2) To calculate the first partial derivative of f with respect to x (var=0) at
point r = [x,y,z], using a FFD with an integration step 1.e-4 at a precision 3**2 use:
fwd_diff_first(*r, var = 0, func = f, delta_x = 1.e-4, order = 3)


RUN

The code can be run directly from command line via:

shell> ./two_body_orbit_int_gen_3d_leap_phys.py <input parameter file>

using the prededfined python interpreter given in the first line,
or:

shell> python3.X two_body_orbit_int_gen_3d_leap_phys.py <input parameter file>

using an alternative python interpreter. Note that a python version 3.X is
required.


INPUT

The input of this program consists of a python (i.e. filename.py) input
parameter file which contains essentially the inital state vectors, i.e.
the initial position r0 and velocities v0 for each body, as well as their
corresponding potential.
The corresponding orbital parameters (e.g. eccentricity, orbital period,
etc.) are then calculated from these.

Note that providing the state vectos r0 and v0 for each body is the preferred
method to specify the initial conditions of the system. However, some methods
are available to impose certain conditions on the orbit, e.g. a circular orbit
*within* the input parameter file (see the example parameter files).


OUTPUT

The output of the code consists of an ascii table redirected to the
directory './output' named after the input parameter file and appended
by the substring "_out.dat".
The table consists of a total of 24 columns. The first 16 contain, in that
order, the time, the specific relative angular momentum, the relative potential
energy, the relative kinetic energy, and (x,vx,y,vt,z,vz) for each of the
bodies. The last 8 columnns contain (x',vx',y',vt',z',vz') for each of the
bodies, where the primed coordinates and velocities correspond to those
on the orbital plane. In other words, these coordinates represent a
rotated version of the intrinsic orbit such that the relative angular momentum
aligns with the z-axis. If the intrinsic orbit has this property,
the primed and unprimed coordinates are identical. Note that the primed z
coordinates are ignored, because they all vanish by definition.
The set of primed coordinates are very useful for checking the numerical result
against the expected analytic solution based on the (potentially osculating)
orbital parameters calculated at runtime.


UNITS

This code is intended mainly for astrophysical applications. Therefore, the
following units are adopted:
  Gravitational constant:   4.301E-06 kpc km^2 / Msun s^2
               Mass unit:       1.000 solMass = 1.988475E+33 g
             Length unit:       1.000 kpc = 3.085678E+21 cm
           Velocity unit:       1.000 km / s = 1.000000E+05 cm / s
               Time unit:       0.978 Gyr = 3.085678E+16 s

Note that the conversion factor of the time is fixed by the others. The
choice of units for G is convenient, since any given potential (or specific
potential energy) automatically has units of (km/s)^2, i.e. identical to the
specific kinetic energy.


VISUALISATION

A simple visualisation of the evolution of the system can be obtained
using the gnuplot script plot_orbit_two_body.gp provided with the code's
distribution.

This script can be run:

1) Directly from the command line via:

shell> gnuplot plot_orbit_two_body.gp

or

2) Within gnuplot via:

gnuplot> load 'plot_orbit_two_body.gp'

The script provides both a full 3D view of the system, or alternatively,
a 2D projection along one of the principal axis of a standard Cartesian
reference frame.

This script accepts a number of input arguments, all of which are set
to reasonable defaults, the obvious exception being the data file.
Other important parameters are the time step, the time output frequency,
and the physical units, all of which are printed to stdout by the
python code at runtime.

Please consult the script's header for additional information on these
input arguments and their default values.


TO DO:


'''

import math
import sys
sys.path.insert(0,"../.")							# include top directory
sys.path.insert(0,"./params")						# include initial conditions directory
import importlib									# needed to import ICs' as module
from ode_int.leapfrog import ode_leap
# from num_diff.forward_diff import fwd_diff_first	# forward finite difference scheme
from num_diff.central_diff import cen_diff_first	# central finite difference scheme (recommended)
from utils import funcs
from astropy import units
import phys_consts as pc

# Collect program argument(s)
if len(sys.argv) < 2:
	print("\nUSAGE:")
	print("{} <input parameter file>\n".format(sys.argv[0]))
	exit()
else:
	ics_file = sys.argv[1]
	ext = ".py"
	if ext in ics_file:
		initialConds = ics_file.replace(ext,'') 	# remove file extension if present
		print("\nRemoved unnecessary extension {} in input parameter file.".format(ext))
	else:
		initialConds = ics_file


# Informative output
print("\nConstants and units:\n")
print("{:>25}{:12.3E} {:15}".format("Gravitational constant:",pc.Grav,"kpc km^2 / Msun s^2"))
massUnit = 1. * units.M_sun
print("{:>25}{:12.3f} = {:6E}".format("Mass unit:",massUnit,massUnit.cgs))
lengthUnit = 1. * units.kpc
print("{:>25}{:12.3f} = {:6E}".format("Length unit:",lengthUnit,lengthUnit.cgs))
velUnit = (1.*units.km/units.s)
print("{:>25}{:12.3f} = {:6E}".format("Velocity unit:",velUnit,velUnit.cgs))
timeUnit = (1.*units.kpc / (units.km / units.s)).decompose().to(units.Gyr)
print("{:>25}{:12.3f} = {:6E}".format("Time unit:",timeUnit,timeUnit.cgs))


# Initial conditions
# Recall: body 1 may be an extended object, while body 2 is a point mass.
# from params import params_ini_3d_two_body_0_bis as ic
ic = importlib.import_module(initialConds)
print("\nGathering input parameters from file {}...".format(initialConds))
t0 = ic.t_0
t1 = ic.t_1
try:
	timeStep = ic.delta_t
except:
	timeStep = 0.001
# Body 1
M1 = ic.Mass1
Phi1 = ic.Potential1
try:
	mass1_cum = ic.Mass1_cum
except:
	print("\nWARNING: No cumulative mass function defined for body 1.")
	print("Assuming a point-like mass distribution.\n")
x10 = ic.x1_0
y10 = ic.y1_0
z10 = ic.z1_0
vx10 = ic.vx1_0
vy10 = ic.vy1_0
vz10 = ic.vz1_0
# Body 2
M2 = ic.Mass2
Phi2 = ic.Potential2
try:
	mass2_cum = ic.Mass2_cum
except:
	print("\nWARNING: No cumulative mass function defined for body 2.")
	print("Assuming a point-like mass distribution.\n")
x20 = ic.x2_0
y20 = ic.y2_0
z20 = ic.z2_0
vx20 = ic.vx2_0
vy20 = ic.vy2_0
vz20 = ic.vz2_0
print("Done.")

# check consistency between potential and cumulative mass function
# TO DO
# if Phi1.__name__ == "Kepler_Pot":

# relative coordinates and velocities
x21_0 = x20 - x10
y21_0 = y20 - y10
z21_0 = z20 - z10
vx21_0 = vx20 - vx10
vy21_0 = vy20 - vy10
vz21_0 = vz20 - vz10
r21_0_vec = [x21_0,y21_0,z21_0]
v21_0_vec = [vx21_0,vy21_0,vz21_0]
r21_0 = funcs.norm(*r21_0_vec)
v21_0 = funcs.norm(*v21_0_vec)

# total mass
Mtot=M1+M2

# Mass of (possibly extended) body 1 at r21_0
try:
	mass1_r21_0 = mass1_cum(*r21_0_vec)
except:
	mass1_r21_0 = M1
# Mass of point-like body 2 at r21_0
try:
	mass2_r21_0 = mass2_cum(*r21_0_vec)
except:
	mass2_r21_0 = M2

# gravitational parameter (only affects the calculation of
# orbital parameters but not the actual orbit calculation)
# grav_param = pc.Grav*Mtot
grav_param = pc.Grav*(mass1_r21_0+mass2_r21_0)

# reduced mass
Mred = (M1*M2)/Mtot											

# Acceleration function definitions (consider moving these into the input parameter file)

# The following functions correspond each to one of the (numerical) partial derivatives of
# either Phi1 or Phi2; the latter must be defined in the input parameter file!

# NOTE: The spatial integration step and the order of accuracy affect
# the conservation of angular momentum, but not the conservation of energy (which
# seems to depend on the time integration scheme only).
#
# Forward finite difference scheme:
# Somewhat counterintuitive, an order 2 integration appears to show a more stable
# evolution of the total angular momentum. Order 3 leads to a gradual increase with
# time, although the error is orders of magnitude smaller than when using order = 2
# for a fixed integration step. Also, the error in the conservation of angular momentum
# is generally orders of magnitude smaller than the error in the conservation of energy (why?).
# Example:
# To calculate the first partial derivative of f with respect to x (var=0) at
# point r = [x,y,z], using an integration step 1.e-4 at a precision 3**2 use:
# fwd_diff_first(*r, var = 0, func = f, delta_x = 1.e-4, order = 3)

# Central finite difference scheme:
# Using an order 2 scheme conserves well both energy and angular momentum.
# An order 4 scheme conserves energy well, and angular momentum to machine precision.
# DO NOT change unless necessary!
intStep = 1.e-4
accOrder = 4

# Note: the array f = (x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
# Recall: Force field = - Grad Phi
def dvx1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_rvec = [_x12,_y12,_z12]
	return  -1.*cen_diff_first(*_rvec, var=0, func=Phi2, delta_x=intStep, order=accOrder)

def dvy1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_rvec = [_x12,_y12,_z12]
	return   -1.*cen_diff_first(*_rvec, var=1, func=Phi2, delta_x=intStep, order=accOrder)

def dvz1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_rvec = [_x12,_y12,_z12]
	return   -1.*cen_diff_first(*_rvec, var=2, func=Phi2, delta_x=intStep, order=accOrder)

def dvx2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_rvec = [_x21,_y21,_z21]
	return  -1.*cen_diff_first(*_rvec, var=0, func=Phi1, delta_x=intStep, order=accOrder)

def dvy2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_rvec = [_x21,_y21,_z21]
	return  -1.*cen_diff_first(*_rvec, var=1, func=Phi1, delta_x=intStep, order=accOrder)

def dvz2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_rvec = [_x21,_y21,_z21]
	return  -1.*cen_diff_first(*_rvec, var=2, func=Phi1, delta_x=intStep, order=accOrder)


# Calculate orbital parameters
# Note that these are intended to characterise the orbit and do not affect the
# orbit integration in any way.

Ltot0_vec = funcs.cross_prod(r21_0_vec,v21_0_vec)
Ltot0 = funcs.norm(*Ltot0_vec)
if Ltot0 > 0:
	semi_latus = funcs.semi_latus_rec(Ltot0,grav_param)
	Ltot0_vec_norm = [ l/Ltot0 for l in Ltot0_vec]
else:
	raise ValueError("Initial conditions imply a purely radial orbit (vanishing angular momentum).")


v_tan_0 = Ltot0 / r21_0
v_rad_0 = math.sqrt(v21_0**2 - v_tan_0**2)

# Specific Laplace-Runge-Lenz vector a.k.a. 'eccentricity vector'
# 'specific' means it is normalised by (G Mtot)*(M_reduced)
ecc_vec = funcs.eccentricity_vec(r21_0_vec,v21_0_vec,Ltot0_vec,grav_param)

# Eccentricity
ecc = funcs.norm(*ecc_vec)

# axes
semimajor_axis = funcs.semimajor(ecc,semi_latus)
semiminor_axis = funcs.semiminor(ecc,semi_latus)
pericen = funcs.pericentre(ecc,semi_latus)
v_peri = Ltot0 / pericen
apocen = funcs.apocentre(ecc,semi_latus)
v_apo = Ltot0 / apocen
orbital_period = funcs.period(semimajor_axis,grav_param)

orb_circum_corr = 0.25 * \
	((semimajor_axis-semiminor_axis)/(semimajor_axis+semiminor_axis))**2
orb_circum_approx = math.pi*(semimajor_axis+semimajor_axis) * \
	(1. + orb_circum_corr + (orb_circum_corr/4.)**4)

orbital_period_peri = orb_circum_approx / v_peri
orbital_period_apo = orb_circum_approx / v_apo

ePot0 = Mred * (Phi1(*r21_0_vec)+Phi2(*r21_0_vec))
eKin0 = Mred * funcs.eKin(*v21_0_vec)

# Determine transformation that maps state vectors onto
# orbital plane; makes use of Rodrigues' rotation formula
# See: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
# The transformation is such that the angular momentum aligns
# with the positive z-axis
z_axis = [0.,0.,1.]
cos_theta = funcs.dot_prod(z_axis,Ltot0_vec_norm)
orbit_incl = math.acos(cos_theta)
k_vec = funcs.cross_prod(Ltot0_vec_norm,z_axis)


# This angle describes to rotation of the LRL vector with respect to the
# x-axis and this determines the orientation of the orbit *in the orbital plane*.
# Need to turn this into a proper function
ecc_vec_rot = funcs.rodrigues_rot(ecc_vec,k_vec,orbit_incl)
apsidal_angle = math.atan2(ecc_vec_rot[1],ecc_vec_rot[0])

print("\nOrbital parameters (osculating for non-Keplerian orbits):\n")
print("{:>40}{:>15}".format("Potential of body 1:",Phi1.__name__))
print("{:>40}{:>15}".format("Potential of body 2:",Phi2.__name__))
print("{:>40}{:15.4E}".format("Total mass of body 1:",M1))
print("{:>40}{:15.4E}".format("Total mass of body 2:",M2))
print("{:>40}{:15.4E}\n".format("Reduced mass:",Mred))

print("{:>40}{:15.4f}".format("Initial rel. distance (r21_0):",r21_0))
print("{:>40}{:15.4f}".format("Initial rel. speed (v21_0):",v21_0))
print("{:>40}{:15.4f}".format("Initial tangential vel. (v_tan_0):",v_tan_0))
print("{:>40}{:15.4f}".format("Initial radial vel. (v_rad_0):",v_rad_0))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("Rel. specific ang. mom. vec. (h_vec):",*Ltot0_vec))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("[normalised]:",*Ltot0_vec_norm))
print("{:>40}{:15.4f}\n".format("Rel. specific ang. mom. (h):",Ltot0))

print("{:>40}{:15.4f}".format("Rel. semi-latus rectum (p):",semi_latus))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("Rel. eccentricity vector (e_vec):",*ecc_vec))
print("{:>40}{:15.4f}".format("Rel. eccentricity (e):",ecc))
print("{:>40}{:15.4f}".format("Rel. semimajor axis (a):",semimajor_axis))
print("{:>40}{:15.4f}".format("Rel. semiminor axis (b):",semiminor_axis))
print("{:>40}{:15.4f}".format("Rel. pericentre (rp):",pericen))
print("{:>40}{:15.4f}".format("Vel. at pericentre (vp):",v_peri))
print("{:>40}{:15.4f}".format("Rel. apocentre (ra):",apocen))
print("{:>40}{:15.4f}".format("Vel. at apocentre (va):",v_apo))
print("{:>40}{:15.4f}".format("Rel. orbital period (T):",orbital_period))
print("{:>40}{:15.4f}".format("Approx. orbit circumference (u):",orb_circum_approx))
print("{:>40}{:15.4f}".format("Approx. pericentric period (Tp):",orbital_period_peri))
print("{:>40}{:15.4f}".format("Approx. apocentric period (Ta):",orbital_period_apo))
print("{:>40}{:15.4f}".format("Apsidal angle (phi_0; deg):",math.degrees(apsidal_angle)))
print("{:>40}{:15.4f}\n".format("Orbital inclination (psi_0; deg):",math.degrees(orbit_incl)))

print("{:>40}{:15.4E}".format("Rel. potential energy (T):",ePot0))
print("{:>40}{:15.4E}".format("Rel. kinetic energy (T):",eKin0))
print("{:>40}{:15.4E}".format("Rel. total energy (T):",ePot0+eKin0))

# Set up time integrator
N = math.ceil((t1 - t0) / timeStep)
ics = [t0, x10, vx10, y10, vy10, z10, vz10, x20, vx20, y20, vy20, z20, vz20]
F = [dvx1dt, dvy1dt, dvz1dt, dvx2dt, dvy2dt, dvz2dt]

print("\nTime range [t0,t1] = [{},{}]".format(t0,t1))
print("Time steps {};  Step size: {:8.2E}".format(N,timeStep))
print("Maximum time step to avoid energy drift: {:12.6E}"\
.format(orbital_period_peri / (math.sqrt(2.)*math.pi)))

# Integrate
# Note: x1 = EoM[0], vx1 = EoM[1], y1 = EoM[2], vy1 = EoM[3], z1 = EoM[4], vz1 = EoM[5],
#       x2 = EoM[6], vx2 = EoM[7], y2 = EoM[8], vy2 = EoM[9], z2 = EoM[10], vz2 = EoM[11]
time, EoM = ode_leap(dr2dt2 = F, rank = 12, initCond = ics, steps = N, stepSize = timeStep)

# Output
outDir = "./output/"

# 1: Intrinsic orbit
outFile = outDir + initialConds + "_out.dat"
outStep = int(N/min(1000,int(1./timeStep)))
print("\nWriting output to file {} with timestep frequency {}\n".format(outFile, outStep))
f = open(outFile, 'wt')
f.write(("{:<9}{:<6}{:8} {:<14}"+"{:6} {:<9}"*2+"{:4} {:<9}"*20+"\n").\
	format("# time ", str(timeUnit.unit), \
		"ang.mom.", str(lengthUnit.unit*velUnit.unit), \
		"ePot", str(velUnit.unit**2), "eKin", str(velUnit.unit**2), \
		"x1", str(lengthUnit.unit), "vx1", str(velUnit.unit), \
		"y1", str(lengthUnit.unit), "vy1", str(velUnit.unit), \
		"z1", str(lengthUnit.unit), "vz1", str(velUnit.unit), \
		"x2", str(lengthUnit.unit), "vx2", str(velUnit.unit), \
		"y2", str(lengthUnit.unit), "vy2", str(velUnit.unit), \
		"z2", str(lengthUnit.unit), "vz2", str(velUnit.unit), \
		"x1_proj", str(lengthUnit.unit), "vx1_proj", str(velUnit.unit), \
		"y1_proj", str(lengthUnit.unit), "vy1_proj", str(velUnit.unit), \
		"x2_proj", str(lengthUnit.unit), "vx2_proj", str(velUnit.unit), \
		"y2_proj", str(lengthUnit.unit), "vy2_proj", str(velUnit.unit)))
eCons = 0.
lCons = 0.
for t in range(0,N,outStep):
	x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2 = \
		EoM[0][t],EoM[1][t],EoM[2][t],EoM[3][t],EoM[4][t],EoM[5][t], \
		EoM[6][t],EoM[7][t],EoM[8][t],EoM[9][t],EoM[10][t],EoM[11][t]
	r1 = [x1,y1,z1]
	v1 = [vx1,vy1,vz1]
	r2 = [x2,y2,z2]
	v2 = [vx2,vy2,vz2]
	r_rel = [x2-x1,y2-y1,z2-z1]
	v_rel = [vx2-vx1,vy2-vy1,vz2-vz1]
# 	Ltot = Mred * funcs.norm(*funcs.cross_prod(r_rel,v_rel))			# rel.ang.mom.
	Ltot = funcs.norm(*funcs.cross_prod(r_rel,v_rel))					# spec.rel.ang.mom.
	ePot1 = Phi2(*r_rel)
	ePot2 = Phi1(*r_rel)
	ePot = Mred * (ePot1+ePot2)											# rel. potential energy (?)
	eKin = Mred * funcs.eKin(*v_rel)									# rel. kin. energy
	r1_rot = funcs.rodrigues_rot(r1,k_vec,orbit_incl)						# rotate vectors
	v1_rot = funcs.rodrigues_rot(v1,k_vec,orbit_incl)
	r2_rot = funcs.rodrigues_rot(r2,k_vec,orbit_incl)
	v2_rot = funcs.rodrigues_rot(v2,k_vec,orbit_incl)
	x1_rot,y1_rot,_ = r1_rot											# ignore z-component
	vx1_rot,vy1_rot,_ = v1_rot
	x2_rot,y2_rot,_ = r2_rot
	vx2_rot,vy2_rot,_ = v2_rot
	f.write(("{:<15.8f}{:<23.10E}"+"{:<16.8E}"*2+"{:<14.4E}"*20+"\n").\
		format(time[t]*timeUnit.value, Ltot, ePot, eKin, \
			x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2, \
			x1_rot, vx1_rot, y1_rot, vy1_rot, x2_rot, vx2_rot, y2_rot, vy2_rot))
	eCons_t = abs((ePot+eKin)/(ePot0+eKin0)-1.)
	if eCons_t > eCons:
		eCons = eCons_t
	lCons_t = abs((Ltot/Ltot0)-1.)
	if lCons_t > lCons:
		lCons = lCons_t
f.close()

print("\n{:>45} {:8.3E} %.".format("Energy conservation to better than",1.e2*eCons))
print("{:>45} {:8.3E} %.\n".format("Angular momentum conservation to better than",1.e2*lCons))

