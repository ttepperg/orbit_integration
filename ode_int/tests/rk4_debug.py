#!/opt/local/bin//python3.5

'''

Author: Thorsten Tepper-Garcia
Date: 02/07/2019

This is an implementation of the Runge-Kutta method of order
4 to integrate a second order ordinary differential equation.

It is intended to become a module.

The background is taken from Kai:

https://math.stackexchange.com/questions/721076/
help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od

'''

def main():

	import math

	def f(x,*y):
		return y[1]

	def g(x,*y):
		return 6*y[0] - y[1]

	def y_ana(x):
		return	math.exp(-3*x) + 2 * math.exp(2*x)

	def dydx_ana(x):
		return	-3*math.exp(-3*x) + 4 * math.exp(2*x)

	def ode_rk4(dydx = None, order = None, initVal = None, steps = None, stepSize = None):

		# Perform sanity checks
		if order is not None:
			k1 = [None] * order
			k2 = [None] * order
			k3 = [None] * order
			k4 = [None] * order
		else:
			raise TypeError("ode_rk4: 'order' is a required argument.")

		if dydx is None:
			raise TypeError("ode_rk4: 'dydx' is a required argument.")

		if len(dydx) != order:
			raise ValueError("ode_rk4: order of 'dydx' does not match 'order'.")

		if steps is None:
			raise TypeError("ode_rk4: 'steps' is a required argument.")

		if stepSize is None:
			raise TypeError("ode_rk4: 'stepSize' is a required argument.")

		# Set initial conditions
		x = [ initVal[0] ]								# independent variable
		y = [ [initVal[m+1]] for m in range(order)]		# m dependent variables		

		print("# i x", end=' ')
		for m in range(order):
			print("rel_err_y{}".format(m), end=' ')
		print()

		# loop over steps
		for s in range(steps):

			# calculate RK coefficients
			Y1 = [y[m][s] for m in range(order)]
			for m in range(order):
				k1[m] = stepSize * dydx[m](x[s], *Y1)

			Y2 = [y[m][s]+0.5*k1[m] for m in range(order)]
			for m in range(order):
				k2[m] = stepSize * dydx[m](x[s]+0.5*stepSize, *Y2)

			Y3 = [y[m][s]+0.5*k2[m] for m in range(order)]
			for m in range(order):
				k3[m] = stepSize * dydx[m](x[s]+0.5*stepSize, *Y3)

			Y4 = [y[m][s]+k3[m] for m in range(order)]
			for m in range(order):
				k4[m] = stepSize * dydx[m](x[s]+stepSize, *Y4)

			# update variables
			x.append(x[s] + stepSize)
			for m in range(order):
				y[m].append(y[m][s] + (1./6.) * (k1[m] + 2*(k2[m]+k3[m]) + k4[m]))

			print("{} {:5.2f}".format(s,x[s]), end=' ')
			for m in range(order):
				print("{:E}".\
				format((y[m][s]-y_ana(x[s]))/y_ana(x[s]),(y[1][s]-dydx_ana(x[s]))/dydx_ana(x[s])), \
				end=' ', flush = True)
			print()


	# Step 0
	a, b = 0., 1.
	NN = 10
	hh = (b-a)/NN

	print("Range [a,b] = [{},{}]".format(a,b))
	print("Steps {};  Step size: {}".format(NN,hh))

	# Step 1
	# Note: the following is the *particular* solution for the
	# given initial conditions:

	# Step 2
	x0 = a
	y0 = 3.
	z0 = 1.

	print("Initial conditions (x0,y0,z0) = ({},{},{})".format(x0,y0,z0))

	# Step 3

	FF = [f,g]
	ic = [x0,y0,z0]

	# function call
	ode_rk4(dydx = FF, order = 2, initVal = ic, steps = NN, stepSize = hh)

if __name__ == "__main__": main()
