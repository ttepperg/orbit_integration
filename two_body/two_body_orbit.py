#!/opt/local/bin//python3.5

'''
Author:	Thorsten Tepper Garcia
Date: 	01/07/2019

See README for information on the code's background, usage, etc.
'''

import math
import sys
sys.path.insert(0,"../.")							# include top directory
sys.path.insert(0,"./init")						# include initial conditions directory
import importlib									# needed to import ICs' as module
from ode_int.leapfrog import ode_leap
# from num_diff.forward_diff import fwd_diff_first	# forward finite difference scheme
from num_diff.central_diff import cen_diff_first	# central finite difference scheme (recommended)
from utils import funcs
# from astropy import units
from config import units
import config.phys_consts as pc


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
print("{:>25}{:12.3f} = {:6E}".format("Mass unit:",units.MASS,units.MASS.cgs))
print("{:>25}{:12.3f} = {:6E}".format("Length unit:",units.LENGTH,units.LENGTH.cgs))
print("{:>25}{:12.3f} = {:6E}".format("Velocity unit:",units.VELOCITY,units.VELOCITY.cgs))
print("{:>25}{:12.3f} = {:6E}".format("Time unit:",units.TIME,units.TIME.cgs))


# Initial conditions
# Recall: body 1 may be an extended object, while body 2 is a point mass.
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


# radial and tangential velocities
v_tan_0 = Ltot0 / r21_0
v_rad_0 = math.sqrt(v21_0**2 - v_tan_0**2)

# specific Laplace-Runge-Lenz vector a.k.a. 'eccentricity vector'
# 'specific' means it is normalised by (G Mtot)*(M_reduced)
ecc_vec = funcs.eccentricity_vec(r21_0_vec,v21_0_vec,Ltot0_vec,grav_param)

# eccentricity
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

# rotate eccentricity vector onto orbital plane
ecc_vec_rot = funcs.rodrigues_rot(ecc_vec,k_vec,orbit_incl)

# The apsidal angle describes to rotation of the eccentricity vector with
# respect to the # x-axis and this determines the orientation of the
# orbit *in the orbital plane*.
apsidal_angle = math.atan2(ecc_vec_rot[1],ecc_vec_rot[0])

# Initia (osculating) orbital parameters of the system
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


# Output to file
outDir = "./output/"
outFile = outDir + initialConds + "_out.dat"
outStep = int(N/min(1000,int(1./timeStep)))
print("\nWriting output to file {} with timestep frequency {}\n".format(outFile, outStep))

f = open(outFile, 'wt')
f.write(("{:<9}{:<6}{:8} {:<14}"+"{:6} {:<9}"*2+"{:4} {:<9}"*20+"\n").\
	format("# time ", str(units.TIME.unit), \
		"ang.mom.", str(units.LENGTH.unit*units.VELOCITY.unit), \
		"ePot", str(units.VELOCITY.unit**2), "eKin", str(units.VELOCITY.unit**2), \
		"x1", str(units.LENGTH.unit), "vx1", str(units.VELOCITY.unit), \
		"y1", str(units.LENGTH.unit), "vy1", str(units.VELOCITY.unit), \
		"z1", str(units.LENGTH.unit), "vz1", str(units.VELOCITY.unit), \
		"x2", str(units.LENGTH.unit), "vx2", str(units.VELOCITY.unit), \
		"y2", str(units.LENGTH.unit), "vy2", str(units.VELOCITY.unit), \
		"z2", str(units.LENGTH.unit), "vz2", str(units.VELOCITY.unit), \
		"x1_proj", str(units.LENGTH.unit), "vx1_proj", str(units.VELOCITY.unit), \
		"y1_proj", str(units.LENGTH.unit), "vy1_proj", str(units.VELOCITY.unit), \
		"x2_proj", str(units.LENGTH.unit), "vx2_proj", str(units.VELOCITY.unit), \
		"y2_proj", str(units.LENGTH.unit), "vy2_proj", str(units.VELOCITY.unit)))
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
	Ltot = funcs.norm(*funcs.cross_prod(r_rel,v_rel))					# spec.rel.ang.mom. (i.e. divided by Mred)
	ePot1 = Phi2(*r_rel)
	ePot2 = Phi1(*r_rel)
	ePot = Mred * (ePot1+ePot2)											# rel. potential energy (?)
	eKin = Mred * funcs.eKin(*v_rel)									# rel. kin. energy
	r1_rot = funcs.rodrigues_rot(r1,k_vec,orbit_incl)					# rotate vectors
	v1_rot = funcs.rodrigues_rot(v1,k_vec,orbit_incl)
	r2_rot = funcs.rodrigues_rot(r2,k_vec,orbit_incl)
	v2_rot = funcs.rodrigues_rot(v2,k_vec,orbit_incl)
	x1_rot,y1_rot,_ = r1_rot											# ignore z-component because is 0.
	vx1_rot,vy1_rot,_ = v1_rot
	x2_rot,y2_rot,_ = r2_rot
	vx2_rot,vy2_rot,_ = v2_rot
	f.write(("{:<15.8f}{:<23.10E}"+"{:<16.8E}"*2+"{:<14.4E}"*20+"\n").\
		format(time[t]*units.TIME.value, Ltot, ePot, eKin, \
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


# experimental
# from utils import io
# io.write_table("test.dat")