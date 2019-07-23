'''
	Author: Thorsten Tepper Garcia

	Original author: Nick Ryan

	Adapted from:

https://nickcdryan.com/2017/09/13/root-finding-algorithms-in-python-line-search-bisection-secant-newton-raphson-boydens-inverse-quadratic-interpolation-brent_root/

	An implementation of Brent's root-finding algorithm.
	
'''

def brent_root(f = None, x0 = None, x1 = None, max_iter=50, tolerance=1e-5):
	"""
	NAME:

		brent_root

	PURPOSE:

		Find a root of the f = 0 using Brent's root-finding algorithm

	INPUT:

		f - a scalar function of x
		x0 - lower bound of searching range
		x1 - upper bound of searching range

	OUTPUT:

		x1 - the root, i.e. f(x1) = 0
          

	HISTORY:
		Unknown - Written - Nick Ryan
		2019-07-20 - Adapted - TTG

	"""
	if f is None:
		raise ValueError("f is required parameter in brent_root")
	elif x0 is None:
		raise ValueError("x0 is required parameter in brent_root")
	elif x1 is None:
		raise ValueError("x1 is required parameter in brent_root")
	else:
		fx0 = f(x0)
		fx1 = f(x1)
 
		assert (fx0 * fx1) <= 0, "Root not bracketed" 
 
		if abs(fx0) < abs(fx1):
			x0, x1 = x1, x0
			fx0, fx1 = fx1, fx0
 
		x2, fx2 = x0, fx0
 
		mflag = True
		steps_taken = 0
 
		while steps_taken < max_iter and abs(x1-x0) > tolerance:
			fx0 = f(x0)
			fx1 = f(x1)
			fx2 = f(x2)
 
			if fx0 != fx2 and fx1 != fx2:
				L0 = (x0 * fx1 * fx2) / ((fx0 - fx1) * (fx0 - fx2))
				L1 = (x1 * fx0 * fx2) / ((fx1 - fx0) * (fx1 - fx2))
				L2 = (x2 * fx1 * fx0) / ((fx2 - fx0) * (fx2 - fx1))
				new = L0 + L1 + L2
 
			else:
				new = x1 - ( (fx1 * (x1 - x0)) / (fx1 - fx0) )
 
			if ((new < ((3 * x0 + x1) / 4) or new > x1) or
				(mflag == True and (abs(new - x1)) >= (abs(x1 - x2) / 2)) or
				(mflag == False and (abs(new - x1)) >= (abs(x2 - d) / 2)) or
				(mflag == True and (abs(x1 - x2)) < tolerance) or
				(mflag == False and (abs(x2 - d)) < tolerance)):
				new = (x0 + x1) / 2
				mflag = True
 
			else:
				mflag = False
 
			fnew = f(new)
			d, x2 = x2, x1
 
			if (fx0 * fnew) < 0:
				x1 = new
			else:
				x0 = new
 
			if abs(fx0) < abs(fx1):
				x0, x1 = x1, x0
 
			steps_taken += 1
 
		return x1



# Test
if __name__ == '__main__':

	f = lambda x: x**3 - 8
 
	root = brent_root(f, 1., 4., tolerance=1.e-5)
	print("root is:", root)
