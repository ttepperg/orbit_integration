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
from class_defs import orbit

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
# (consider moving this into body class)
ic = importlib.import_module(initialConds)
print("\nGathering input parameters from file {}...".format(initialConds))

try:
	t0 = ic.t_0
except:
	t0 = 0.
try:
	t1 = ic.t_1
except:
	t1 = 10.
try:
	timeStep = ic.delta_t
except:
	timeStep = 0.001

# Body 1
try:
	Pot = ic.Potential1
except:
	Pot = None
r10_vec = [ic.x1_0,ic.y1_0,ic.z1_0]
v10_vec = [ic.vx1_0,ic.vy1_0,ic.vz1_0]
body1 = body.Body(mass=ic.Mass1,pot=Pot,r_vec=r10_vec,v_vec=v10_vec)

# Body 2
try:
	Pot = ic.Potential2
except:
	Pot = None
r20_vec = [ic.x2_0,ic.y2_0,ic.z2_0]
v20_vec = [ic.vx2_0,ic.vy2_0,ic.vz2_0]
body2 = body.Body(mass=ic.Mass2,pot=Pot,r_vec=r20_vec,v_vec=v20_vec)

print("Done.")

# Initialise an Orbit object
orbit = orbit.Orbit(body1,body2)

# Output to stdout
# (consider redirecting these to output file as well)
orbit.orbit_info()

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
ics = \
[t0, body1.x, body1.vx, body1.y, body1.vy, body1.z, body1.vz, body2.x, body2.vx, body2.y, body2.vy, body2.z, body2.vz]
F = [dvx1dt, dvy1dt, dvz1dt, dvx2dt, dvy2dt, dvz2dt]

print("\nTime range [t0,t1] = [{},{}]".format(t0,t1))
print("Time steps {};  Step size: {:8.2E}".format(N,timeStep))
print("Maximum time step to avoid energy drift: {:12.6E}".format(orbit.energy_drift_lim()))

# Integrate orbit
# Note: x1 = EoM[0], vx1 = EoM[1], y1 = EoM[2], vy1 = EoM[3], z1 = EoM[4], vz1 = EoM[5],
#       x2 = EoM[6], vx2 = EoM[7], y2 = EoM[8], vy2 = EoM[9], z2 = EoM[10], vz2 = EoM[11]
# 		each EoM[i] is 'function' of time, i.e. an array where each item corresponds to a
# 		given integration time, i.e. EoM[i][t]
time, EoM = ode_leap(dr2dt2 = F, rank = 12, initCond = ics, steps = N, stepSize = timeStep)


# initial relative energies
ePot0 = orbit.energy_pot()
eKin0 = orbit.energy_kin()

# initial relative angular momentum
Ltot0 = orbit.ang_mom()

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
	ePot = orbit.mred() * (ePot1+ePot2)											# rel. potential energy (?)
	eKin = orbit.mred() * funcs.eKin(*v_rel)									# rel. kin. energy
	r1_rot = orbit.orbital_plane(r1)									# rotate vectors onto orbital plane
	v1_rot = orbit.orbital_plane(v1)
	r2_rot = orbit.orbital_plane(r2)
	v2_rot = orbit.orbital_plane(v2)
	x1_rot,y1_rot,_ = r1_rot											# ignore z-component because is 0.
	vx1_rot,vy1_rot,_ = v1_rot
	x2_rot,y2_rot,_ = r2_rot
	vx2_rot,vy2_rot,_ = v2_rot

	f.write(("{:<15.8f}{:<23.10E}"+"{:<16.8E}"*2+"{:<14.4E}"*22+"\n").\
		format(time[t]*units.TIME.value, Ltot, ePot, eKin, \
			x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2, \
			x1_rot, vx1_rot, y1_rot, vy1_rot, x2_rot, vx2_rot, y2_rot, vy2_rot, \
			*orbit.kepler_orbit_cartesian(t*2.*math.pi/N)))
			
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