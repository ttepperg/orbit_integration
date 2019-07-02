#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019

Integrate the second order ordinary differential equation:

d2y/dx2 + dy/dx - 6 y  = 0

with solution:

y(x) = exp[-3x] + 2 exp[2x]

for initial conditions y(0) = 3 and dy/dx(0) = 1.

For an thorough explanation see ode_int/test_rk4.py

'''

import math
import sys
sys.path.insert(0,'../../.')
from ode_int.rk4 import ode_rk4

# Function definitions
# Transform 2nd order ODE as a system of two 1st order ODE
def f1(x,*y):
	return y[1]

def f2(x,*y):
	return 6*y[0] - y[1]

# Analytic solution
# Note: the following is the *particular* solution for the
# given initial conditions:
def sol(x):
	return	math.exp(-3*x) + 2 * math.exp(2*x)


# Initial conditions
f1_0 = 3.
f2_0 = 1.

# Define integration range, steps, and step size
a, b = 0., 1.
N = 10
h = (b-a)/N

print("Range [a,b] = [{},{}]".format(a,b))
print("Steps {};  Step size: {}".format(N,h))

# Set up integrator
x_0 = a
ic = (x_0,f1_0,f2_0)
F = [f1,f2]

print("Initial conditions (x0,y0,z0) = {}".format(ic))

# ODE integration
# Note: the actual solution will be stored y_num[0], and its derivative in y_num[1]
t, y_num = ode_rk4(dydx = F, order = 2, initVal = ic, steps = N, stepSize = h)
y_ana = sol(1.)

print("y(1)_numeric = {}".format(y_num[0][N]))
print("y(1)_analytic = {}".format(y_ana))
print("Absolute relative error: {}".format(abs(y_num[0][N] - y_ana)/y_ana))
