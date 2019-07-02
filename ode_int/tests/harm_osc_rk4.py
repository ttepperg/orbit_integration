#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019

Integrate the second order ordinary differential equation
for the harmonic oscillator:

d2y/dt2 + w^2 y  = 0

with solution:

y(x) = A * sin(wt+phi)

for initial conditions y(0) = A sin(phi) and dy/dx(0) = A w cos(phi).

using the 4th order Runge-Kutta method.

'''

import math
import sys
sys.path.insert(0,'../../.')
from ode_int.rk4 import ode_rk4

# Frequency
omega = 3.1

# Initial conditions
y_0 = 1.2
dydx_0 = 0.4

# Define integration range, steps, and step size
a, b = 0., 2*math.pi
N = 100
h = (b-a)/N

print("Range [a,b] = [{},{}]".format(a,b))
print("Steps {};  Step size: {}".format(N,h))

# Function definitions
# Transform 2nd order ODE as a system of two 1st order ODE
# f1 = dydt = z
# f2 = dzdt = -w^2 y
# Note: the array f = (y,z), i.e y = f[0], z = f[1]
def f1(t,*f):
	return f[1]

def f2(t,*f):
	return -omega**2 * f[0]

# Analytic solution
def sol(t,a,w,p):
	return	a * math.sin(w * t + p)

# Set up integrator
t_0 = a
ic = (t_0,y_0,dydx_0)
F = [f1,f2]

print("Initial conditions (x0,y0,z0) = {}".format(ic))

# ODE integration
# Note: the actual solution will be stored y_num[0], and its derivative in y_num[1]
time, y_num = ode_rk4(dydx = F, order = 2, initVal = ic, steps = N, stepSize = h)

# Parameters of the analytic solution
if dydx_0 != 0:
	phase = math.atan(omega * y_0/dydx_0)
else:
	phase = 0.5*math.pi
if math.sin(phase) != 0:
	Amp = y_0 / math.sin(phase)
else:
	Amp = dydx_0 / omega

##############################################################################
# Write to ASCII file

outDir = './output/'
fileOut = outDir+"harm_osc_rk4.dat"

# write table
print("Writing output to ascii file {}".format(fileOut))
f = open(fileOut, 'wt')
f.write("# Amplitude (A): {}\n".format(Amp))
f.write("# Phase (p;rad): {}\n".format(phase))
f.write("{:<12} {:<12} {:<12}\n".format("# t","y_num(x)","A sin(wt+p)"))
for t in range(N):
	f.write("{:<12.6f} {:<12.6f} {:<12.6f}\n".\
	format(time[t],y_num[0][t], sol(time[t],Amp,omega,phase)))
f.close()

##############################################################################
print("Done.\n")

exit()