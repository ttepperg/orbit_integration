#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019

This is a test implementation of the Runge-Kutta method of order
4 to integrate a second order ordinary differential equation:

d2y/dx2 + dy/dx - 6 y  = 0

with solution:

y(x) = exp[-3x] + 2 exp[2x]

for initial conditions y(0) = 3 and dy/dx(0) = 1.

The example and solution are taken from:

https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od

'''

# BACKGROUND

# Step 0: Chose a range for the solution, a number of steps, and the
# step size
# Range: a < x < b
# Steps: N
# Step size: h = (b - a) / N

# Step 1: Transform the 2nd order ODE into a system of coupled
# 1st order ODEs:
# f = dy/dx = z
# g = dz/dx = 6 y - z


# Step 2: Define the initial conditions x0, y0 = y(x0) and z0 = z(y0).

# Step 3: Apply the RK4 method, for a 2nd order ODE.
# NOTE 1: Because this is a system of ODEs, the RK4 steps need to be
# performed in a particular order (need reference!).
# NOTE 2: It is relatively straightforward to generalise the
# method to a Nth order ODE, provided it can be transformed
# into a system of N possibly coupled 1st order ODEs.
# See link above.

# k1 = h * f(x0,y0,z0)
# l1 = h * g(x0,y0,z0)
# k2 = h * f(x0+0.5*h, y0+0.5*k1, z0+0.5*l1)
# l2 = h * g(x0+0.5*h, y0+0.5*k1, z0+0.5*l1)
# k3 = h * f(x0+0.5*h, y0+0.5*k2, z0+0.5*l2)
# l3 = h * g(x0+0.5*h, y0+0.5*k2, z0+0.5*l2)
# k4 = h * f(x0+h, y0+k3, z0+l3)
# l4 = h * g(x0+h, y0+k3, z0+l3)
# y1 = y0 + (1./6.) * (k1 + 2*(k2+k3) + k4)
# z1 = z0 + (1./6.) * (l1 + 2*(l2+l3) + l4)

# IMPLEMENTATION

import math

# Step 0
a, b = 0., 1.
N = 10
h = (b-a)/N

print("Range [a,b] = [{},{}]".format(a,b))
print("Steps {};  Step size: {}".format(N,h))

# Step 1
def f(x,y,z):
	return z

def g(x,y,z):
	return 6*y - z

# Note: the following is the *particular* solution for the
# given initial conditions:
def y_ana(x):
	return	math.exp(-3*x) + 2 * math.exp(2*x)

def dydx_ana(x):
	return	-3*math.exp(-3*x) + 4 * math.exp(2*x)

# Step 2
x0 = a
y0 = 3.
z0 = 1.

print("Initial conditions (x0,y0,z0) = ({},{},{})".format(x0,y0,z0))

# Step 3
x = []
y = []
z = []

x.append(x0)
y.append(y0)
z.append(z0)

print("# i x rel_err_y rel_err_z")
for i in range(N+1):				# recall that range() is exclusive
	k1 = h * f(x[i],y[i],z[i])
	l1 = h * g(x[i],y[i],z[i])
	k2 = h * f(x[i]+0.5*h, y[i]+0.5*k1, z[i]+0.5*l1)
	l2 = h * g(x[i]+0.5*h, y[i]+0.5*k1, z[i]+0.5*l1)
	k3 = h * f(x[i]+0.5*h, y[i]+0.5*k2, z[i]+0.5*l2)
	l3 = h * g(x[i]+0.5*h, y[i]+0.5*k2, z[i]+0.5*l2)
	k4 = h * f(x[i]+h, y[i]+k3, z[i]+l3)
	l4 = h * g(x[i]+h, y[i]+k3, z[i]+l3)
	x.append(x[i] + h)
	y.append(y[i] + (1./6.) * (k1 + 2*(k2+k3) + k4))
	z.append(z[i] + (1./6.) * (l1 + 2*(l2+l3) + l4))
	print("{} {:5.2f} {:E} {:E}".format(i,x[i],(y[i]-y_ana(x[i]))/y_ana(x[i]),(z[i]-dydx_ana(x[i]))/dydx_ana(x[i])))
