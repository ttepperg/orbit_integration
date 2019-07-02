'''
(C) Thor Tepper Garcia 2019

This module contains an implementation of the Runge-Kutta method of order
4 with constant step size to integrate an possibly vectorial, m-order ordinary
differential equation (ODE), provided it can be expressed as a system of m 1st order
ODEs.

The background is taken from Kai:

https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od


Last modified: 18/06/2019

DOCUMENTATION:

To Do


'''

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

	if len(initVal) != order+1:
		raise \
		ValueError("ode_rk4: expected order + 1 initial conditions, got {}'.".format(len(initVal)))

	if steps is None:
		raise TypeError("ode_rk4: 'steps' is a required argument.")

	if stepSize is None:
		raise TypeError("ode_rk4: 'stepSize' is a required argument.")

	# Set initial conditions
	x = [ initVal[0] ]								# independent variable
	y = [ [initVal[m+1]] for m in range(order)]		# m dependent variables		

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

	return x, y
