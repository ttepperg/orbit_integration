'''
(C) Thor Tepper Garcia 2019


Implementation of a central finite difference scheme to calculate
the first (partial) derivative of an arbitrary function to accuracy
of order 2 or 4.

The background is taken from:


https://en.wikipedia.org/wiki/Finite_difference_coefficient

'''

def cen_diff_first(*r0, var = None, func = None, delta_x = None, order = None):
	"""
		Returns the first (partial) derivative of f = func with
		respect to its var variable r0 to given order of accuracy in
		the integration step delta_x, i.e. ~O(delta_X**order).
		The parameter var takes values from 0 to len(r0) - 1.
		Currently only accuracies of order 2 or 3 are implemented.
	"""
	funcName = "cen_diff_first"

	# sanity check
	if len(r0) <= 0:
		raise TypeError("expected r0 of length at least one 1 in {}, got {}.".format(funcName,len(r0)))
	if var is None:
		raise TypeError("var is a required parameter in {}.".format(funcName))
	if var > len(r0)-1:
		raise TypeError("var ({}) exceeds the largest index of r0 ({}) in {}.".format(var,len(r0)-1,funcName))
	if func is None:
		raise TypeError("func is a required parameter in {}.".format(funcName))
	if delta_x is None:
		raise TypeError("delta_x is a required parameter in {}.".format(funcName))
	if order is None:
		raise TypeError("order is a required parameter in {}.".format(funcName))

	# set finite difference coefficients and stencil depending on order
	if order == 2:
		coeff = [-0.5, 0.5]
		r_1 = list(r0)
		r1 = list(r0)
		r_1[var] = r0[var] - delta_x
		r1[var] = r0[var] + delta_x

		# error check
		try: func(*r_1)
		except:
			raise ValueError("func not defined at r_1={} in {}.".format(r_1,funcName))
		try: func(*r1)
		except:
			raise ValueError("func not defined at r1={} in {}.".format(r1,funcName))

		if delta_x > 0.:
			delta_x_inv = 1./delta_x
			return delta_x_inv*(coeff[0]*func(*r_1)+coeff[1]*func(*r1))
		else:
			raise ValueError("Non-valid delta_x value ({}) in {}.".format(delta_x,funcName))

	elif order == 4:
		coeff = [1./12., -2./3., 2./3., -1./12.]
		r_2 = list(r0)
		r_1 = list(r0)
		r1 = list(r0)
		r2 = list(r0)
		r_2[var] = r0[var] - 2 * delta_x
		r_1[var] = r0[var] - delta_x
		r1[var] = r0[var] + delta_x
		r2[var] = r0[var] + 2 * delta_x

		# error check
		try: func(*r_2)
		except:
			raise ValueError("func not defined at r_2={} in {}.".format(r_2,funcName))
		try: func(*r_1)
		except:
			raise ValueError("func not defined at r_1={} in {}.".format(r_1,funcName))
		try: func(*r1)
		except:
			raise ValueError("func not defined at r1={} in {}.".format(r1,funcName))
		try: func(*r2)
		except:
			raise ValueError("func not defined at r2={} in {}.".format(r2,funcName))

		if delta_x > 0.:
			delta_x_inv = 1./delta_x
			return delta_x_inv*(coeff[0]*func(*r_2)+coeff[1]*func(*r_1)+coeff[2]*func(*r1)+coeff[3]*func(*r2))
		else:
			raise ValueError("Non-valid delta_x value ({}) in {}.".format(delta_x,funcName))

	else:
		raise ValueError("currently order can only be 2 or 3, got {} in {}.".format(order, funcName))


	