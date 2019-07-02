#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019

The 3D equations of motion of two body of mass m orbiting around each
other potential are:

dr1dt = v1
dv1dt = - Grad Phi2

dr2dt = v2
dv2dt = - Grad Phi1

Here,

r := position vector of body 1 or 2
v := velocity vector of body 1 or 2
Grad := gradient operator (nabla)
Phi(r) = := potential of body 1 or 2

Here, we consider the bodies to be represented by point masses, such
that Phi is given by a Kepler potential of a point mass M at distance |r|:

Phi(r) = - G M / |r| .

dx1dt = vx1
dy1dt = vy1
dz1dt = vz1
dvx1dt = - (G M2 / |r12|^2) ((x1-x2)/|r12|)
dvy1dt = - (G M2 / |r12|^2) ((y1-y2)/|r12|)
dvz1dt = - (G M2 / |r12|^2) ((z1-z2)/|r12|)

dx2dt = vx2
dy2dt = vy2
dz2dt = vz2
dvx2dt = - (G M1 / |r12|^2) ((x2-x1)/|r12|)
dvy2dt = - (G M1 / |r12|^2) ((y2-y1)/|r12|)
dvz2dt = - (G M1 / |r12|^2) ((z2-z1)/|r12|)

Here, r12 = r1 - r2. Note that the force is radial and directed towards
the centre of the potential at all times.


Even though the equation of motion is second order in dt,
due to its 3D vector nature, it is equivalent to a system of six
first order ODEs for each body. Thus, the rank of the system is determined
by the order m and the dimensionality of the system D as 2 * m * D = 12 in
this case.

Here, the system of equations is numerically integrated using a 2nd order
leapfrog (kick-drift-kick) method. Note that this does conserve energy well,
even over significant integration times, but may lead to orbit precession it
the integration step is not small enough.
Angular momentum appears to be conserved exactly, i.e. to machine precision.

'''

import math
import sys
sys.path.insert(0,"../.")	# include top directory
from ode_int.leapfrog import ode_leap
from utils import funcs
from astropy import units

import phys_consts as pc


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

# Function definitions

# Acceleration
# Note: the array f = (x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
# Important: force along line connecting the bodies, directed towards the
# opposite body.
# NOTE: these functions correspond to each of the partial derivatives of
# Phi1 or Phi2; the latter are defined in the input parameter file (see below).
def dvx1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_r2 = funcs.norm2(_x12,_y12,_z12)
	_r = funcs.norm(_x12,_y12,_z12)
	return  -1. * (pc.Grav * M2 / _r2) * (_x12/_r)

def dvy1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_r2 = funcs.norm2(_x12,_y12,_z12)
	_r = funcs.norm(_x12,_y12,_z12)
	return  -1. * (pc.Grav * M2 / _r2) * (_y12/_r)

def dvz1dt(t,*f):
	_x12 = f[0]-f[6]
	_y12 = f[2]-f[8]
	_z12 = f[4]-f[10]
	_r2 = funcs.norm2(_x12,_y12,_z12)
	_r = funcs.norm(_x12,_y12,_z12)
	return  -1. * (pc.Grav * M2 / _r2) * (_z12/_r)

def dvx2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_r2 = funcs.norm2(_x21,_y21,_z21)
	_r = funcs.norm(_x21,_y21,_z21)
	return  -1. * (pc.Grav * M1 / _r2) * (_x21/_r)

def dvy2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_r2 = funcs.norm2(_x21,_y21,_z21)
	_r = funcs.norm(_x21,_y21,_z21)
	return  -1. * (pc.Grav * M1 / _r2) * (_y21/_r)

def dvz2dt(t,*f):
	_x21 = f[6]-f[0]
	_y21 = f[8]-f[2]
	_z21 = f[10]-f[4]
	_r2 = funcs.norm2(_x21,_y21,_z21)
	_r = funcs.norm(_x21,_y21,_z21)
	return  -1. * (pc.Grav * M1 / _r2) * (_z21/_r)


# Initial conditions
from params import params_ini_3d_two_body_2 as ic
t0 = ic.t_0
t1 = ic.t_1
M1 = ic.Mass1
Phi1 = ic.Potential1
x10 = ic.x1_0
y10 = ic.y1_0
z10 = ic.z1_0
vx10 = ic.vx1_0
vy10 = ic.vy1_0
vz10 = ic.vz1_0
M2 = ic.Mass2
Phi2 = ic.Potential2
x20 = ic.x2_0
y20 = ic.y2_0
z20 = ic.z2_0
vx20 = ic.vx2_0
vy20 = ic.vy2_0
vz20 = ic.vz2_0

# sanity check
if Phi1.__name__ != "Kepler_Pot" or Phi1.__name__ != "Kepler_Pot":
	raise TypeError("Phi1 and Phi2 can only be Kepler_Potential")

# Calculate orbital parameters
x210 = x20 - x10
y210 = y20 - y10
z210 = z20 - z10
vx210 = vx20 - vx10
vy210 = vy20 - vy10
vz210 = vz20 - vz10
r210 = funcs.norm(x210,y210,z210)
v210 = funcs.norm(vx210,vy210,vz210)
semimajor_axis = funcs.semimajor_from_vis(r210,v210**2,amp = pc.Grav*(M1+M2))

ang_mom_spec = funcs.cross_prod(x210,vx210,y210,vy210,z210,vz210)
radius_p = funcs.mean_radius(ang_mom_spec, amp = pc.Grav*(M1+M2))
ecc = funcs.eccentricity(radius_p,semimajor_axis)
semiminor_axis = funcs.semiminor(radius_p,semimajor_axis)
orbital_period = funcs.period(semimajor_axis, amp = pc.Grav*(M1+M2))
pericen = funcs.pericentre(semimajor_axis,ecc)
apocen = funcs.apocentre(semimajor_axis,ecc)

print("\nOrbital parameters")
print("{:>30}{:12.4E}".format("Mass of body 1:",M1))
print("{:>30}{:12.4E}".format("Mass of body 2:",M2))
print("{:>30}{:12.4f}".format("Initial rel. distance (r210):",r210))
print("{:>30}{:12.4f}".format("Initial rel. speed (v210):",v210))
print("{:>30}{:12.4f}".format("Rel. eccentricity (e):",ecc))
print("{:>30}{:12.4f}".format("Rel. semimajor axis (a):",semimajor_axis))
print("{:>30}{:12.4f}".format("Rel. semiminor axis (b):",semiminor_axis))
print("{:>30}{:12.4f}".format("Rel. pericentre (rp):",pericen))
print("{:>30}{:12.4f}".format("Rel. apocentre (ra):",apocen))
print("{:>30}{:12.4f}".format("Rel. radius (p):",radius_p))
print("{:>30}{:12.4f}".format("Rel. specific ang. mom. (h):",ang_mom_spec))
print("{:>30}{:12.4f}".format("Rel. orbital period (T):",orbital_period))


# Set up integrator
timeStep = 0.001
N = math.ceil((t1 - t0) / timeStep)
ics = [t0, x10, vx10, y10, vy10, z10, vz10, x20, vx20, y20, vy20, z20, vz20]
F = [dvx1dt, dvy1dt, dvz1dt, dvx2dt, dvy2dt, dvz2dt]

print("\nRange [t0,t1] = [{},{}]".format(t0,t1))
print("Steps {};  Step size: {}".format(N,timeStep))

# Integrate
# Note: x1 = EoM[0], vx1 = EoM[1], y1 = EoM[2], vy1 = EoM[3], z1 = EoM[4], vz1 = EoM[5],
#       x2 = EoM[6], vx2 = EoM[7], y2 = EoM[8], vy2 = EoM[9], z2 = EoM[10], vz2 = EoM[11]
time, EoM = ode_leap(dr2dt2 = F, rank = 12, initCond = ics, steps = N, stepSize = timeStep)

# Output
outDir = "./output/"
outFile = outDir + "two_body_orbit_int_kepler_3d_leap_phys.dat"
outStep = int(N/1000)
print("\nWriting output to file {} with timestep frequency {}\n".format(outFile, outStep))
f = open(outFile, 'wt')
f.write(("{:<9}{:<6}{:8} {:<14}" + "{:4} {:<9}"*14 + "\n").\
	format("# time ", str(timeUnit.unit), \
		"ang.mom.", str(lengthUnit.unit*velUnit.unit), \
		"epot", str(velUnit.unit**2), "ekin", str(velUnit.unit**2), \
		"x1", str(lengthUnit.unit), "vx1", str(velUnit.unit), \
		"y1", str(lengthUnit.unit), "vy1", str(velUnit.unit), \
		"z1", str(lengthUnit.unit), "vz1", str(velUnit.unit), \
		"x2", str(lengthUnit.unit), "vx2", str(velUnit.unit), \
		"y2", str(lengthUnit.unit), "vy2", str(velUnit.unit), \
		"z2", str(lengthUnit.unit), "vz2", str(velUnit.unit)))
Mtot = M1+M2
Mred = M1*M2/Mtot
for t in range(0,N,outStep):
	x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2 = \
		EoM[0][t],EoM[1][t],EoM[2][t],EoM[3][t],EoM[4][t],EoM[5][t], \
		EoM[6][t],EoM[7][t],EoM[8][t],EoM[9][t],EoM[10][t],EoM[11][t]
	Rcom_x, Rcom_y, Rcom_z = \
		(M1*x1+M2*x2)/Mtot, (M1*y1+M2*y2)/Mtot, (M1*z1+M2*z2)/Mtot			# position of centre of mass
	Vcom_x, Vcom_y, Vcom_z = \
		(M1*vx1+M2*vx2)/Mtot, (M1*vy1+M2*vy2)/Mtot, (M1*vz1+M2*vz2)/Mtot	# velocity of centre of mass
	Lcom = \
		Mtot*funcs.cross_prod(Rcom_x,Vcom_x,Rcom_y,Vcom_y,Rcom_z,Vcom_z)	# magnitude of ang.mom. of CoM
	L12 = \
		Mred*funcs.cross_prod(x1-x2,vx1-vx2,y1-y2,vy1-vy2,z1-z2,vz1-vz2)	# magnitude of ang.mom. wrt CoM
	V1 = M1*funcs.ePot(x1-x2, y1-y2, z1-z2, pot = Phi2)						# potential energy of body 1
	V2 = M2*funcs.ePot(x1-x2, y1-y2, z1-z2, pot = Phi1)						# potential energy of body 2
	Tcom = Mtot*funcs.eKin(Vcom_x,Vcom_y,Vcom_z)							# kinetic energy of CoM
	T12 = Mred*funcs.eKin(vx1-vx2,vy1-vy2,vz1-vz2)							# kinetic energy wrt CoM
	f.write(("{:<15.8f}{:<23.10E}"+"{:<14.4E}"*14+"\n").\
		format(time[t]*timeUnit.value, Lcom+L12, 0.5*(V1+V2), Tcom+T12, \
			x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2)) #(*)
f.close()

# (*) Note: The total potential due to internal forces is 0.5 the sum of potential pairs (see Goldstein, eq. 1.36)
