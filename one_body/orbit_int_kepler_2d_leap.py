#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019

The 2D equation of motion of a body of mass m orbiting in a potential
Phi is:

drdt = v
dvdt = - Grad Phi

Here,

r := position vector
v := velocity vector
Grad := gradient operator (nabla)
Phi(r) = - G M / |r|

i.e. a Kepler potential of a point mass M at distance |r|.

Note that the force is radial and directed towards the centre of the
potential at all times.

In components, with r = (x,y) and v=(vx,vy)

dxdt = vx
dydt = vy
dvxdt = - (G M / |r|^2) (x/|r|)
dvydt = - (G M / |r|^2) (y/|r|)


I.e. even though the equation of motion is second order in dt,
due to its vector nature (2D), it is equivalent to a system of four
first order ODEs. Thus, the rank of the system is determined by the
order m and the dimensionality of the system D as m * D.

Here, the system of equations is numerically integrated using a 2nd order
leapfrog (kick-drift-kick) method. Note that this does conserve energy well,
even over significant integration times, but leads to orbit precession it
the integration step is not small.

For an identical code but using physical units, see orbit_int_kepler_2d_phys.py.

For an alternative integration method, see orbit_int_kepler_2d_rk4.py.

NOTE: See doc/exercise_sheet.pdf for a meaningful set of initial conditions.

'''

import math
import sys
sys.path.insert(0,"../.")	# include the top directory
from ode_int.leapfrog import ode_leap

# Gravitational constant (in natural units for now):
Grav = 1.

# Mass (units?)
M = 10.

# Function definitions

def norm2(x,y):
	return (x**2 + y**2)

def norm(x,y):
	return math.sqrt(norm2(x,y))

# Note: the array f = (x,vx,y,vy)
def dvxdt(t,*f):
	_r2 = norm2(f[0],f[2])
	_r = norm(f[0],f[2])
	return  -1. * (Grav * M / _r2) * (f[0]/_r)

def dvydt(t,*f):
	_r2 = norm2(f[0],f[2])
	_r = norm(f[0],f[2])
	return  -1. * (Grav * M / _r2) * (f[2]/_r)

def Phi(*r):
	_r = norm(*r)
	return  -1. * (Grav * M / _r)

def Epot(*r):
	return Phi(*r)

def Ekin(*v):
	return 0.5 * norm2(*v)


# Initial conditions (units?)
t_0 = 0.
x_0 = 1.
y_0 = 0.
dxdt_0 = 0.
dydt_0 = 1.

# Set up integrator
t_1 = 10
h = 0.001
N = math.ceil((t_1 - t_0) / h)
ics = [t_0, x_0, dxdt_0, y_0, dydt_0]
F = [dvxdt, dvydt]

print("Range [t0,t1] = [{},{}]".format(t_0,t_1))
print("Steps {};  Step size: {}".format(N,h))

# Integrate
# Note: x = EoM[0], vx = EoM[1], y = EoM[2], vy = EoM[3]
time, EoM = ode_leap(dr2dt2 = F, rank = 4, initCond = ics, steps = N, stepSize = h)

# Output
outDir = "./output/"
outFile = outDir + "orbit_int_kepler_2d_leap.dat"
outStep = int(N/1000)
f = open(outFile, 'wt')
f.write("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".\
format("# time","ekin","epot","x(t)","vx(t)","y(t)","vy(t)"))
for t in range(0,N,10):
	epot = Epot(EoM[0][t], EoM[2][t])
	ekin = Ekin(EoM[1][t], EoM[3][t])
	f.write("{:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f} {:<12.6f}\n".\
	format(time[t], epot, ekin, EoM[0][t], EoM[1][t], EoM[2][t], EoM[3][t]))
f.close()
