#!/opt/local/bin//python3.5

'''
Author:			Thorsten Tepper Garcia
Last modified:	04/07/2019

See README for information on the code's background, usage, etc.

'''

# Import necessary modules

# built-in
import math
import importlib									# needed to import ICs' as module
import sys

# costum modules
sys.path.insert(0,"../.")							# include top directory in module search path
sys.path.insert(0,"./init")							# include initial conditions directory
from ode_int.leapfrog import ode_leap
# from num_diff.forward_diff import fwd_diff_first	# forward finite difference scheme
from num_diff.central_diff import cen_diff_first	# central finite difference scheme (recommended)
from utils import funcs
from utils import io
from config import units
import config.phys_consts as pc
from class_defs import body

# Collect input argument(s)
initialConds = io.get_input()

# Informative output
# (consider moving this into the units.py module)
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
r10_vec = [ic.x1_0,ic.y1_0,ic.z1_0]
v10_vec = [ic.vx1_0,ic.vy1_0,ic.vz1_0]
body1 = body.Body(mass=ic.Mass1,pot=ic.Potential1,r_vec=r10_vec,v_vec=v10_vec)

# Body 2
r20_vec = [ic.x2_0,ic.y2_0,ic.z2_0]
v20_vec = [ic.vx2_0,ic.vy2_0,ic.vz2_0]
body2 = body.Body(mass=ic.Mass2,pot=ic.Potential2,r_vec=r20_vec,v_vec=v20_vec)

print("Done.")

# relative coordinates and velocities
r21_0_vec = body2.pos_rel(body1)
v21_0_vec = body2.vel_rel(body1)
r21_0 = body2.dist_rel(body1)
v21_0 = body2.speed_rel(body1)


# gravitational parameter (only affects the calculation of
# orbital parameters but not the actual orbit calculation)
# uses masses of bodies at r21_0_vec
grav_param = pc.Grav*(body1.mass_cum(*r21_0_vec)+body2.mass_cum(*r21_0_vec))

# total mass (NOT necessarily equal to body1.mass_cum(*r21_0_vec)+body2.mass_cum(*r21_0_vec)!)
Mtot=body1.mass+body2.mass

# reduced mass
Mred = (body1.mass*body2.mass)/Mtot											

# Calculate orbital parameters
# NOTE: These are intended to characterise the orbit and do not affect the
# orbit integration in any way. They are, however, very useful when it comes
# to a comparison between the code's result and the analytic solution, as well
# and to check the conservation laws (total energy and angular momentum).
# (consider encapsulating all of this in a class)

Ltot0_vec = funcs.cross_prod(r21_0_vec,v21_0_vec)
Ltot0 = funcs.norm(*Ltot0_vec)
if Ltot0 > 0:
	Ltot0_vec_norm = [ l/Ltot0 for l in Ltot0_vec]
else:
	raise ValueError("Initial conditions imply a purely radial orbit (vanishing angular momentum).")


# radial and tangential velocities
v_tan_0 = Ltot0 / r21_0
v_rad_0 = math.sqrt(v21_0**2-v_tan_0**2)

# specific Laplace-Runge-Lenz vector a.k.a. 'eccentricity vector'
# 'specific' means it is normalised by (G Mtot)*(M_reduced)
ecc_vec = funcs.eccentricity_vec(r21_0_vec,v21_0_vec,Ltot0_vec,grav_param)

# eccentricity
ecc = funcs.norm(*ecc_vec)

# principal axes
semi_latus = funcs.semi_latus_rec(Ltot0,grav_param)
semimajor_axis = funcs.semimajor(ecc,semi_latus)
semiminor_axis = funcs.semiminor(ecc,semi_latus)
pericen = funcs.pericentre(ecc,semi_latus)
apocen = funcs.apocentre(ecc,semi_latus)

v_peri = Ltot0 / pericen
v_apo = Ltot0 / apocen
orbital_circum = funcs.circumference(semimajor_axis,semiminor_axis)
orbital_period = funcs.period(semimajor_axis,grav_param)
orbital_period_peri = orbital_circum / v_peri
orbital_period_apo = orbital_circum / v_apo

ePot0 = Mred * (body1.potential(*r21_0_vec)+body2.potential(*r21_0_vec))
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

# rotate angular momentum vector onto z-axis
Ltot0_vec_rot = funcs.rodrigues_rot(Ltot0_vec,k_vec,orbit_incl)

# rotate eccentricity vector onto orbital plane
ecc_vec_rot = funcs.rodrigues_rot(ecc_vec,k_vec,orbit_incl)

# The apsidal angle describes to rotation of the eccentricity vector,
# which points from apocentre to pericentre, with respect to the x-axis
# it determines the orientation of the orbit *in the orbital plane*
apsidal_angle = math.atan2(ecc_vec_rot[1],ecc_vec_rot[0])

# Output to stdout (consider including these in output file)
print("\nInitial (osculating) orbital parameters of the system:\n")
print("{:>40}{:>15}".format("Potential of body 1:",body1.potential.__name__))
print("{:>40}{:>15}".format("Potential of body 2:",body2.potential.__name__))
print("{:>40}{:15.4E}".format("Total mass of body 1:",body1.mass))
print("{:>40}{:15.4E}".format("Total mass of body 2:",body2.mass))
print("{:>40}{:15.4E}\n".format("Reduced mass:",Mred))

print("{:>40}{:15.4f}".format("Initial rel. distance (r21_0):",r21_0))
print("{:>40}{:15.4f}".format("Initial rel. speed (v21_0):",v21_0))
print("{:>40}{:15.4f}".format("Initial tangential vel. (v_tan_0):",v_tan_0))
print("{:>40}{:15.4f}".format("Initial radial vel. (v_rad_0):",v_rad_0))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("Rel. specific ang. mom. vec. (h_vec):",*Ltot0_vec))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("[normalised]:",*Ltot0_vec_norm))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("[along z-axis]:",*Ltot0_vec_rot))
print("{:>40}{:15.4f}\n".format("Rel. specific ang. mom. (h):",Ltot0))

print("{:>40}{:15.4f}".format("Rel. semi-latus rectum (p):",semi_latus))
print("{:>40} ({:5.3f},{:5.3f},{:5.3f})".format("Rel. eccentricity vector (e_vec):",*ecc_vec))
print("{:>40}{:15.4f}".format("Rel. eccentricity (e):",ecc))
print("{:>40}{:15.4f}".format("Rel. semimajor axis (a):",semimajor_axis))
print("{:>40}{:15.4f}".format("Rel. semiminor axis (b):",semiminor_axis))
print("{:>40}{:15.4f}".format("Apsidal angle (phi_0; deg):",math.degrees(apsidal_angle)))
print("{:>40}{:15.4f}".format("Orbital inclination (psi_0; deg):",math.degrees(orbit_incl)))
print("{:>40}{:15.4f}".format("Rel. pericentre (rp):",pericen))
print("{:>40}{:15.4f}".format("Vel. at pericentre (vp):",v_peri))
if ecc < 1.:
	print("{:>40}{:15.4f}".format("Rel. apocentre (ra):",apocen))
	print("{:>40}{:15.4f}".format("Vel. at apocentre (va):",v_apo))
	print("{:>40}{:15.4f}".format("Rel. orbital period (T):",orbital_period))
	print("{:>40}{:15.4f}".format("Approx. orbit circumference (u):",orbital_circum))
	print("{:>40}{:15.4f}".format("Approx. pericentric period (Tp):",orbital_period_peri))
	print("{:>40}{:15.4f}\n".format("Approx. apocentric period (Ta):",orbital_period_apo))

print("{:>40}{:15.4E}".format("Rel. potential energy (T):",ePot0))
print("{:>40}{:15.4E}".format("Rel. kinetic energy (T):",eKin0))
print("{:>40}{:15.4E}".format("Rel. total energy (T):",ePot0+eKin0))



# Define acceleration function definitions
# (consider moving these into the input parameter file or a function or a class)

# The following functions correspond each to one of the (numerical) partial derivatives of
# either body1.potential or body2.potential; the latter must be defined in the input parameter file!

# Central finite difference scheme parameters:
# 2nd order scheme conserves well both energy and angular momentum.
# 4th order scheme conserves energy well, and angular momentum generally to
# machine precision.
# DO NOT change the values of the following parameters unless absolutely necessary!
intStep = 1.e-4
accOrder = 4

# Note: the array f = (x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
# Recall: Force field = - Grad Phi
def dvx1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_rvec = [_x12,_y12,_z12]
	return  -1.*cen_diff_first(*_rvec, var=0, func=body2.potential, delta_x=intStep, order=accOrder)

def dvy1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_rvec = [_x12,_y12,_z12]
	return   -1.*cen_diff_first(*_rvec, var=1, func=body2.potential, delta_x=intStep, order=accOrder)

def dvz1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_rvec = [_x12,_y12,_z12]
	return   -1.*cen_diff_first(*_rvec, var=2, func=body2.potential, delta_x=intStep, order=accOrder)

def dvx2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_rvec = [_x21,_y21,_z21]
	return  -1.*cen_diff_first(*_rvec, var=0, func=body1.potential, delta_x=intStep, order=accOrder)

def dvy2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_rvec = [_x21,_y21,_z21]
	return  -1.*cen_diff_first(*_rvec, var=1, func=body1.potential, delta_x=intStep, order=accOrder)

def dvz2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_rvec = [_x21,_y21,_z21]
	return  -1.*cen_diff_first(*_rvec, var=2, func=body1.potential, delta_x=intStep, order=accOrder)

# Set up time integrator
N = math.ceil((t1 - t0) / timeStep)
ics = [t0, body1.x, body1.vx, body1.y, body1.vy, body1.z, body1.vz, body2.x, body2.vx, body2.y, body2.vy, body2.z, body2.vz]
F = [dvx1dt, dvy1dt, dvz1dt, dvx2dt, dvy2dt, dvz2dt]

print("\nTime range [t0,t1] = [{},{}]".format(t0,t1))
print("Time steps {};  Step size: {:8.2E}".format(N,timeStep))
if not math.isnan(orbital_period_peri):
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

f.write(("{:<9}{:<6}{:8} {:<14}"+"{:6} {:<9}"*2+"{:4} {:<9}"*22+"\n").\
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
		"y2_proj", str(units.LENGTH.unit), "vy2_proj", str(units.VELOCITY.unit), \
		"KeplerOrbitAna_X", str(units.LENGTH.unit), "KeplerOrbitAna_Y", str(units.LENGTH.unit)))

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
	ePot1 = body2.potential(*r_rel)
	ePot2 = body1.potential(*r_rel)
	ePot = Mred * (ePot1+ePot2)											# rel. potential energy (?)
	eKin = Mred * funcs.eKin(*v_rel)									# rel. kin. energy
	r1_rot = funcs.rodrigues_rot(r1,k_vec,orbit_incl)					# rotate vectors onto orbital plane
	v1_rot = funcs.rodrigues_rot(v1,k_vec,orbit_incl)
	r2_rot = funcs.rodrigues_rot(r2,k_vec,orbit_incl)
	v2_rot = funcs.rodrigues_rot(v2,k_vec,orbit_incl)
	x1_rot,y1_rot,_ = r1_rot											# ignore z-component because is 0.
	vx1_rot,vy1_rot,_ = v1_rot
	x2_rot,y2_rot,_ = r2_rot
	vx2_rot,vy2_rot,_ = v2_rot

	f.write(("{:<15.8f}{:<23.10E}"+"{:<16.8E}"*2+"{:<14.4E}"*22+"\n").\
		format(time[t]*units.TIME.value, Ltot, ePot, eKin, \
			x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2, \
			x1_rot, vx1_rot, y1_rot, vy1_rot, x2_rot, vx2_rot, y2_rot, vy2_rot, \
			*funcs.kepler_orbit_cartesian(t*2.*math.pi/N,semi_latus,ecc,apsidal_angle)))

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
# io.write_table("test.dat")