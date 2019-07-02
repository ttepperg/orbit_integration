#!/opt/local/bin//python3.5

'''
(C) Thor Tepper Garcia 2019


Test the calculation of a functions (partial) derivative using a central finite
difference scheme of order 2 or 4.

'''

import sys
sys.path.insert(0,'../../.')
from utils import funcs
from central_diff import cen_diff_first


def main():

# Test cases

	f = funcs.Kepler_Potential(amp = 1.)

	# analytic derivative of Kepler_Potential in 1D
	def fdiff_ana_1d(x = None, amp = None):
		if x > 0 and amp is not None:
			return amp / x**2
		else:
			raise ValueError("Invalid parameters x or amp in fdiff_ana_1d.")


	# analytic partial derivative of Kepler_Potential in 2D wrt x
	def fpartdiff_ana_x(*rvec, amp = None):
		_x = rvec[0]
		_y = rvec[1]
		_r2 = funcs.norm2(_x,_y)
		_r = funcs.norm(_x,_y)
		return  (amp / _r2) * (_x/_r)

	# analytic partial derivative of Kepler_Potential in 2D wrt y
	def fpartdiff_ana_y(*rvec, amp = None):
		_x = rvec[0]
		_y = rvec[1]
		_r2 = funcs.norm2(_x,_y)
		_r = funcs.norm(_x,_y)
		return  (amp / _r2) * (_y/_r)


	# order of accuracy (2 or 3)
	m = 4
	print("Num.Deriv.   Ana.Deriv.   Abs.Rel.Error  Exp.Acc.")
	for i in range(1,20):
		x = i*0.01
		der_num = cen_diff_first(x, var = 0, func = f, delta_x = 0.0001, order = m)
		der_ana = fdiff_ana_1d(x, 1.)
		print("{:15.4E}   {:15.4E}   {:15.4E}  {:15.4E}".format(der_num,der_ana,abs((der_num/der_ana)-1.), 0.001**m))
	print()

	print("Part.Num.Deriv.x   Part.Ana.Deriv.x     Abs.Rel.Error      Exp. Acc.")
	for i in range(1,20):
		x = i*0.01
		y = i*1.2
		r = [x,y]
		der_num = cen_diff_first(*r, var = 0, func = f, delta_x = 0.001, order = m)
		der_ana = fpartdiff_ana_x(*r, amp = 1.)
		print("{:15.4E}   {:15.4E}   {:15.4E}  {:15.4E}".format(der_num,der_ana,abs((der_num/der_ana)-1.), 0.001**m))
	print()

	print("Part.Num.Deriv.y   Part.Ana.Deriv.y      Abs.Rel.Error      Exp. Acc.")
	for i in range(1,20):
		x = -i*0.01
		y = i*1.2
		r = [x,y]
		der_num = cen_diff_first(*r, var = 1, func = f, delta_x = 0.001, order = m)
		der_ana = fpartdiff_ana_y(*r, amp = 1.)
		print("{:15.4E}   {:15.4E}   {:15.4E}  {:15.4E}".format(der_num,der_ana,abs((der_num/der_ana)-1.), 0.001**m))
	print()


if __name__ == '__main__': main()
