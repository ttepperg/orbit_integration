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
from utils import funcs
from utils import io
from class_defs import body
from class_defs import orbit
from config import units


# Collect input argument
initialConds = io.get_input()

# Set initial conditions
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

# Output initial orbital parameters to stdout
# (consider redirecting these to output file as well)
orbit.orbit_info()

# Integrate orbit
time_steps, time, EoM = orbit.integrate(t0,t1,timeStep)

# Initial relative energies
ePot0 = orbit.energy_pot()
eKin0 = orbit.energy_kin()

# Initial relative angular momentum
Ltot0 = orbit.ang_mom()

# Output to file
outDir = "./output/"
outFile = outDir + initialConds + "_out.dat"
outStep = int(time_steps/min(1000,int(1./timeStep)))

f = open(outFile, 'wt')

f.write(("{:<9}{:<6}{:8} {:<14}"+"{:6} {:<9}"*2+"{:4} {:<9}"*22+"\n").\
	format("# time ", units.TIME_UNIT_STR, \
		"ang.mom.", units.ANG_MOM_UNIT_STR, \
		"ePot", units.ENERGY_UNIT_STR, "eKin", units.ENERGY_UNIT_STR, \
		"x1", units.LENGTH_UNIT_STR, "vx1", units.VEL_UNIT_STR, \
		"y1", units.LENGTH_UNIT_STR, "vy1", units.VEL_UNIT_STR, \
		"z1", units.LENGTH_UNIT_STR, "vz1", units.VEL_UNIT_STR, \
		"x2", units.LENGTH_UNIT_STR, "vx2", units.VEL_UNIT_STR, \
		"y2", units.LENGTH_UNIT_STR, "vy2", units.VEL_UNIT_STR, \
		"z2", units.LENGTH_UNIT_STR, "vz2", units.VEL_UNIT_STR, \
		"x1_proj", units.LENGTH_UNIT_STR, "vx1_proj", units.VEL_UNIT_STR, \
		"y1_proj", units.LENGTH_UNIT_STR, "vy1_proj", units.VEL_UNIT_STR, \
		"x2_proj", units.LENGTH_UNIT_STR, "vx2_proj", units.VEL_UNIT_STR, \
		"y2_proj", units.LENGTH_UNIT_STR, "vy2_proj", units.VEL_UNIT_STR, \
		"KeplerOrbitAna_X", units.LENGTH_UNIT_STR, "KeplerOrbitAna_Y", units.LENGTH_UNIT_STR))

eCons = 0.
lCons = 0.

for t in range(0,time_steps,outStep):

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
			*orbit.kepler_orbit_cartesian(t*2.*math.pi/time_steps)))
			
	eCons_t = abs((ePot+eKin)/(ePot0+eKin0)-1.)
	if eCons_t > eCons:
		eCons = eCons_t
	lCons_t = abs((Ltot/Ltot0)-1.)
	if lCons_t > lCons:
		lCons = lCons_t

f.close()

print("\n{:>45} {:8.3E} %.".format("Energy conservation to better than",1.e2*eCons))
print("{:>45} {:8.3E} %.\n".format("Angular momentum conservation to better than",1.e2*lCons))


print("\nOutput written to file {} with timestep frequency {}\n".format(outFile, outStep))
