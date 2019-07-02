#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019


The 3D equation of motion of a body of mass m orbiting in the potential
Phi of a central mass M is:

drdt = v
dvdt = - Grad Phi

Here,

r := position vector
v := velocity vector
Grad := gradient operator (nabla)
Phi(r) := potential

Here, we consider point masses, such that Phi is given by
a Kepler potential of a point mass M at distance |r|:

Phi(r) = - G M / |r| .

Note that the force is radial and directed towards the centre of the
potential at all times.

In components, with r = (x,y) and v=(vx,vy)

dxdt = vx
dydt = vy
dydt = vz
dvxdt = - (G M / |r|^2) (x/|r|)
dvydt = - (G M / |r|^2) (y/|r|)
dvzdt = - (G M / |r|^2) (z/|r|)


I.e. even though the equation of motion is second order in dt,
due to its 3D vector nature, it is equivalent to a system of six
first order ODEs. Thus, the rank of the system is determined by the
order m and the dimensionality of the system D as m * D.

Here, the system of equations is numerically integrated using a 2nd order
leapfrog (kick-drift-kick) method. Note that this does conserve energy well,
even over significant integration times, but may lead to orbit precession it
the integration step is not small enough.
Angular momentum appears to be conserved exactly, i.e. to machine precision.

NOTE: See doc/exercise_sheet.pdf for a meaningful set of initial conditions.

'''

import math
import sys
sys.path.insert(0,"../.")	# include the top directory
from ode_int.leapfrog import ode_leap
from utils import funcs
from astropy import units

import config.phys_consts as pc


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
# Note: the array f = (x,vx,y,vy,z,vz)
def dvxdt(t,*f):
	_r2 = funcs.norm2(*[w for w in f[::2]])
	_r = funcs.norm(*[w for w in f[::2]])
	return  -1. * (pc.Grav * M / _r2) * (f[0]/_r)

def dvydt(t,*f):
	_r2 = funcs.norm2(*[w for w in f[::2]])
	_r = funcs.norm(*[w for w in f[::2]])
	return  -1. * (pc.Grav * M / _r2) * (f[2]/_r)

def dvzdt(t,*f):
	_r2 = funcs.norm2(*[w for w in f[::2]])
	_r = funcs.norm(*[w for w in f[::2]])
	return  -1. * (pc.Grav * M / _r2) * (f[4]/_r)


# Initial conditions
from params import params_ini_3d as ic
t0 = ic.t_0
t1 = ic.t_1
M = ic.Mass
Phi = ic.Potential
x0 = ic.x_0
y0 = ic.y_0
z0 = ic.z_0
vx0 = ic.vx_0
vy0 = ic.vy_0
vz0 = ic.vz_0

# Orbital parameters
r0_vec = [x0,y0,z0]
v0_vec = [vx0,vy0,vz0]
r0 = funcs.norm(*r0_vec)
v0 = funcs.norm(*v0_vec)

try:
	semimajor_axis
except:
	semimajor_axis = funcs.semimajor_from_vis(r0,v0,amp = pc.Grav*M)

ang_mom_spec = funcs.norm(*funcs.cross_prod(r0_vec,v0_vec))
radius_p = funcs.semi_latus_rec(ang_mom_spec, amp = pc.Grav*M)
ecc = math.sqrt(1. - radius_p/semimajor_axis)
semiminor_axis = funcs.semiminor(radius_p,semimajor_axis)
orbital_period = funcs.period(semimajor_axis, amp = pc.Grav*M)

print("\nOrbital parameters")
print("{:>25}{:12.4f}".format("Initial distance (r0):",r0))
print("{:>25}{:12.4f}".format("Initial speed (v0):",v0))
print("{:>25}{:12.4f}".format("Eccentricity (e):",ecc))
print("{:>25}{:12.4f}".format("Semimajor axis (a):",semimajor_axis))
print("{:>25}{:12.4f}".format("Semiminor axis (b):",semiminor_axis))
print("{:>25}{:12.4f}".format("Radius (p):",radius_p))
print("{:>25}{:12.4f}".format("Specific ang. mom. (h):",ang_mom_spec))
print("{:>25}{:12.4f}".format("Orbital period (T):",orbital_period))


# Set up integrator
h = 0.001
N = math.ceil((t1 - t0) / h)
ics = [t0, x0, vx0, y0, vy0, z0, vz0]
F = [dvxdt, dvydt, dvzdt]

print("\nRange [t0,t1] = [{},{}]".format(t0,t1))
print("Steps {};  Step size: {}".format(N,h))

# Integrate
# Note: x = EoM[0], vx = EoM[1], y = EoM[2], vy = EoM[3], z = EoM[4], vz = EoM[5]
time, EoM = ode_leap(dr2dt2 = F, rank = 6, initCond = ics, steps = N, stepSize = h)


# Output
outDir = "./output/"
outFile = outDir + "orbit_int_kepler_3d_leap_phys.dat"
outStep = int(N/1000)
print("\nWriting output to file {} with timestep frequency {}\n".format(outFile, outStep))
f = open(outFile, 'wt')
f.write("{:<9}{:<6}{:4} {:<9}{:4} {:<9}{:8} {:<14}{:4} {:<9}{:4} {:<9}{:4} {:<9}{:4} {:<9}{:4} {:<9}{:4} {:<9}\n".\
format("# time ", str(timeUnit.unit), \
	"epot", str(velUnit.unit**2), "ekin", str(velUnit.unit**2), \
	"ang.mom.", str(lengthUnit.unit*velUnit.unit), \
	"x", str(lengthUnit.unit), "vy", str(velUnit.unit), \
	"y", str(lengthUnit.unit), "vx", str(velUnit.unit), \
	"z", str(lengthUnit.unit), "vz", str(velUnit.unit)))
for t in range(0,N,outStep):
	r = [EoM[0][t], EoM[2][t], EoM[4][t]]
	vel = [EoM[1][t], EoM[3][t], EoM[5][t]]
	V = funcs.ePot(*r, pot = Phi)
	T = funcs.eKin(*vel)
	L = funcs.norm(*funcs.cross_prod(r,vel))
	f.write("{:<15.6f}{:<14.6f}{:<14.6f}{:<23.10E}{:<14.6f}{:<14.6f}{:<14.6f}{:<14.6f}{:<14.6f}{:<14.6f}\n".\
	format(time[t]*timeUnit.value, V, T, L, \
		EoM[0][t], EoM[1][t], EoM[2][t], EoM[3][t], EoM[4][t], EoM[5][t]))
f.close()

