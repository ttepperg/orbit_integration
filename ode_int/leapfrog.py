'''
(C) Thor Tepper Garcia 2019

This module contains an implementation of the leapfrog method (kick-drift-kick)
with constant step size to integrate an possibly vectorial (dimension D), mth order ordinary
differential equation (ODE), provided it can be expressed as a system of m * D 1st order
ODEs. The rank of the system is hence m * D.

The background is taken from:

https://en.wikipedia.org/wiki/Leapfrog_integration

Last modified: 18/06/2019

DOCUMENTATION:

To Do

'''

def ode_leap(dr2dt2 = None, rank = None, initCond = None, steps = None, stepSize = None):

	# Perform sanity checks
	if rank is not None:
		rank2 = int(rank/2)
	else:
		raise TypeError("ode_leap: 'rank' is a required argument.")

	if (rank % 2) != 0:
		raise TypeError("ode_rk4: 'rank' is is expected to be an even number.")

	if dr2dt2 is None:
		raise TypeError("ode_leap: 'dr2dt2' is a required argument.")

	if len(dr2dt2) != rank2:
		raise ValueError("ode_leap: rank of 'dr2dt2' ({}) does not match 'rank/2' ({}) .".format(len(dr2dt2),rank2))

	if len(initCond) != rank+1:
		raise \
		ValueError("ode_leap: expected 'rank + 1' = {} initial conditions, got {}'.".format(rank+1,len(initCond)))

	if steps is None:
		raise TypeError("ode_leap: 'steps' is a required argument.")

	if stepSize is None:
		raise TypeError("ode_leap: 'stepSize' is a required argument.")

	# Set initial conditions
	t = [ initCond[0] ]								# independent variable
	y = [ [initCond[m+1]] for m in range(rank)]		# m dependent variables		

	# loop over steps
	for s in range(steps):

		# The following assumes that:
		# - y[0], y[2], ... contain the coordinates
		# - y[1], y[3], ... contain their corresponding velocity(ies)
		# - dr2dt2[0] to dr2dt2[rank/2] contain the acceleration(s)

		# kick: update velocities at first half-step
		Y = [y[m][s] for m in range(rank)]
		for m in range(1,rank,2):
			n = int((m-1)/2)
			y[m].append( y[m][s] + 0.5*stepSize * dr2dt2[n](t[s],*Y) )

		# drift: update coordinates at full step 
		for m in (range(0,rank,2)):
			y[m].append( y[m][s] + stepSize * y[m+1][s+1] )

		# kick: update velocities at second half-step
		Y = [y[m][s+1] for m in range(rank)]
		for m in range(1,rank,2):
			n = int((m-1)/2)
			y[m][s+1] = y[m][s+1] + 0.5*stepSize * dr2dt2[n](t[s],*Y)

		# update independent variable
		t.append(t[s] + stepSize)

	return t, y
