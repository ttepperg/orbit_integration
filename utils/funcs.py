'''
Author:	Thorsten Tepper Garcia
Date:	26/06/2019
'''

import math
import sys
from config.phys_consts import Grav, Infinity
from num_diff.central_diff import cen_diff_first
sys.path.insert(0,'../../.')
from ode_int.rk4 import ode_rk4
from roots.brent import brent_root

# METRICS AND TRANSFORMATIONS

def norm2(*r):
	'''Euclidean distance squared'''
	return sum(w**2 for w in r)


def norm(*r):
	'''Euclidean distance'''
	return math.sqrt(norm2(*r))


def cross_prod(r = None, v = None):
	'''Cross product between vectors r and v, i.e. r x v = - v x r'''
	if r is None:
		raise ValueError("r is a required parameter in cross_prod")
	elif v is None:
		raise ValueError("v is a required parameter in cross_prod")
	elif len(r) != len(v):
		raise ValueError("r and v must have the same length in cross_prod")
	elif len(r) < 2:
		raise ValueError("cannot calculate the cross product in less than 2D in cross_prod")
	else:
		lx = 0.
		ly = 0.
		if len(r) > 2:
			lx = r[1]*v[2]-r[2]*v[1]
			ly = -r[0]*v[2]+r[2]*v[0]
		lz = r[0]*v[1]-r[1]*v[0]
		return lx,ly,lz


def dot_prod(r = None, v = None):
	'''Dot product between vectors r and v, i.e. r . v = v . r'''
	if r is None:
		raise ValueError("r is a required parameter in dot_prod")
	elif v is None:
		raise ValueError("v is a required parameter in dot_prod")
	elif len(r) != len(v):
		raise ValueError("r and v must have the same length in dot_prod")
	elif len(r) < 1:
		raise ValueError("cannot calculate the dot product in less than 1D in dot_prod")
	else:
		dp = 0.
		for i in range(len(r)):
			dp += r[i]*v[i]
		return dp


def rodrigues_rot(vec = None, axis = None, incl = None):
	'''Implements Rodrigues' rotation formula to rotate a vector 'vec'
	by an angle 'incl' around a prescribed 'axis' represented by an unit
	vector k:

	vec_rot = vec cos_theta + (k x vec) sin_theta + k(k . vec)(1-cos_theta)

	See: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

	'''
	if vec is None:
		raise ValueError("vec is a required parameter in rodrigues_rot")
	elif axis is None:
		raise ValueError("axis is a required parameter in rodrigues_rot")
	elif incl is None:
		raise ValueError("incl is a required parameter in rodrigues_rot")
	elif len(vec) != len(axis):
		raise ValueError("vec and axis must have the same length in rodrigues_rot")
	else:
		# ensure axis vector is normalised
		k_norm = norm(*axis)
		if k_norm > 0.:
			k_vec = [k / k_norm for k in axis]
			kcrossv = cross_prod(k_vec,vec)
			kdotv = dot_prod(k_vec,vec)
			cos_theta = math.cos(incl)
			sin_theta = math.sin(incl)
			vec_rot = \
				[a*cos_theta + b*sin_theta + c*kdotv*(1.-cos_theta) \
					for a,b,c in zip(vec,kcrossv,k_vec)]
			return vec_rot
		else:
			# vec and axis are parallel
			return vec


def grad(*r, func = None):
	'''Returns the gradient of a scalar function at r=[x,y,z].
	Grad func is calculated using a central finite difference
	scheme of order 4 and constant step of size 1.e-4.'''
	if func is None:
		raise ValueError("func is a required parameter in grad")
	elif len(r) < 1:
		raise ValueError("r must have dimensionality of at least 1 in grad, but is {}".format(len(r)))
	else:
		grad_func = []
		for i in range(len(r)):
			grad_func.append(cen_diff_first(*r, var=i, func=func, delta_x=1.e-4, order=4))
		return grad_func


def grad_r(*r, func = None):
	'''Returns the directional derivative of a scalar function along u=[x,y,z]/|r|
	at r=[x,y,z], given by the scalar producto between u and Grad func.
	u is the unitary vector along r.
	Grad func is calculated using a central finite difference
	scheme of order 4 and constant step of size 1.e-4.'''
	if func is None:
		raise ValueError("func is a required parameter in grad")
	elif len(r) < 1:
		raise ValueError("r must have dimensionality of at least 1 in grad, but is {}".format(len(r)))
	else:
		grad_func = grad(*r,func=func)
		_r = norm(*r)
		if _r <= 0:
			raise ValueError("r must be have a positive magnitude in grad_r")
		else:
			u = [x/_r for x in r]
			return dot_prod(u,grad_func)



# POTENTIALS, MASS PROFILES, ETC.

# KEPLER MODEL

# Kepler potential
def Kepler_Potential(mass = None):
	'''Implements a wrapper for the Kepler potential where amp = G*M.'''
	if mass is None:
		raise ValueError("mass is a required argument in Kepler_Potential")
	elif mass <= 0:
		raise ValueError("mass must be positive in Kepler_Potential")
	else:
		def Kepler_Pot(*r):
			_amp = Grav*mass
			_r = norm(*r)
			if _r > 0:
				return -1. * (_amp / _r)
			else:
				raise ValueError("Zero or negative radius in Kepler_Potential")
		return Kepler_Pot


# Kepler mass
# Note that this function is trivial, but is included here
# and used in the code for completeness and self-consistency.
def Kepler_Mass(mass = None):
	if mass is None:
		raise ValueError("mass is a required parameter in Kepler_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in Kepler_Mass")
	else:
		def Kepler_M(*r):
			return mass
		return Kepler_M



# HOMOGENEOUS SPHERE

# Homogeneous sphere potential
def Sphere_Potential(mass = None, a = None):
	'''Implements a wrapper for the potential of an homogeneous sphere of given mass and radius a. Reduces to Kepler potential for r >= a'''
	if mass is None:
		raise ValueError("mass is a required argument in Sphere_Potential")
	elif mass <= 0:
		raise ValueError("mass must be positive in Sphere_Potential")
	elif a is None:
		raise ValueError("a is a required argument in Sphere_Potential")
	elif a < 0:
		raise ValueError("a must be non-negative in Sphere_Potential")
	else:
		def Sphere_Pot(*r):
			_amp = Grav * mass
			_r = norm(*r)
			_r2 = _r**2
			_a2 = a**2
			if _r < a:
				return 0.5 * _amp / (a * _a2) * (_r2 - 3. * _a2)
			else:
				return -1 * _amp / _r

		# private attributes (= parent func. params) to allow for access from outside
		Sphere_Pot._mass = mass
		Sphere_Pot._a = a
		return Sphere_Pot


# Homogeneous sphere density profile
def Sphere_Density(mass = None, a = None):
	'''Returns the Homogeneous sphere density at r, which is constant for r < a and zero otherwise'''
	if mass is None:
		raise ValueError("mass is a required parameter in Sphere_Density")
	elif mass <= 0:
		raise ValueError("mass must be positive in Sphere_Density")
	elif a is None:
		raise ValueError("a is a required parameter in Sphere_Density")
	elif a <= 0:
		raise ValueError("a must be positive in Sphere_Density")
	else:
		def Sphere_dens(*r):
			_r = norm(*r)
			if _r < a:
				return 3. * mass / (4. * math.pi * a**3)
			else:
				return 0.
		return Sphere_dens


# Homogeneous sphere cumulative mass
def Sphere_Mass(mass = None, a = None):
	'''Returns the cumulative mass of an homogeneous sphere at r'''
	if mass is None:
		raise ValueError("mass is a required parameter in Sphere_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in Sphere_Mass")
	elif a is None:
		raise ValueError("a is a required parameter in Sphere_Mass")
	elif a < 0:
		raise ValueError("a must be non-negative in Sphere_Mass")
	else:
		def Sphere_M(*r):
			_amp = mass / a**3
			_r = norm(*r)
			if _r < a:
				return _amp * _r**3
			else:
				return mass
		return Sphere_M


# Homogeneous sphere velocity dispersion
def Sphere_VelDisp(mass = None, a = None):
	'''Returns the 1D Homogeneous sphere radial velocity dispersion at r, sigma_r(r)^2. Can be easily derived from the spherically symmetric Jeans equation for isotropic systems (beta = 0) and constant density. In this case, the velocity dispersion is equal to minus the potential.
	IMPORTANT: Strictly speaking, the vel.disp. should be 0 for r > a, but I won't implement that explicity as it may cause a division by 0 when calculating the dynamical friction force. Instead, I set its value to a very large value so as to nullify the dyn.fric. force
	'''
	if mass is None:
		raise ValueError("mass is a required argument in Sphere_VelDisp")
	elif mass <= 0:
		raise ValueError("mass must be positive in Sphere_VelDisp")
	elif a is None:
		raise ValueError("a is a required argument in Sphere_VelDisp")
	elif a < 0:
		raise ValueError("a must be non-negative in Sphere_VelDisp")
	else:
		def Sphere_veldisp(*r):
			_r = norm(*r)
			if _r < a:
				return -1 * Sphere_Potential(mass,a)
			else:
				return 1.e100 # should be 0, but this yield a 0 dyn.fric. force

		return Sphere_veldisp



# PLUMMER (1911) MODEL

# Plummer potential
def Plummer_Potential(mass = None, a = None):
	'''Implements a wrapper for the Plummer (1911) potential.
	Reduces to Kepler potential for a=0.'''
	if mass is None:
		raise ValueError("mass is a required argument in Plummer_Potential")
	elif mass <= 0:
		raise ValueError("mass must be positive in Plummer_Potential")
	elif a is None:
		raise ValueError("a is a required argument in Plummer_Potential")
	elif a < 0:
		raise ValueError("a must be non-negative in Plummer_Potential")
	else:
		def Plummer_Pot(*r):
			_amp = Grav*mass
			_r2 = norm2(*r)
			_d = math.sqrt(_r2+a**2)
			if _d > 0:
				return -1. * (_amp / _d)
			else:
				raise ValueError("Zero or negative denominator in Plummer_Potential")
		# private attributes (= parent func. params) to allow for access from outside
		Plummer_Pot._mass = mass
		Plummer_Pot._a = a
		return Plummer_Pot


# Plummer density profile
def Plummer_Density(mass = None, a = None):
	'''Returns the Plummer density at r'''
	if mass is None:
		raise ValueError("mass is a required parameter in Plummer_Density")
	elif mass <= 0:
		raise ValueError("mass must be positive in Plummer_Density")
	elif a is None:
		raise ValueError("a is a required parameter in Plummer_Density")
	elif a <= 0:
		raise ValueError("a must be positive in Plummer_Density")
	else:
		def Plummer_dens(*r):
			_amp = 3.*mass*a**2 / (4.*math.pi)
			_r = norm(*r)
			_d2 = (a**2 + _r**2)
			if _d2 > 0:
				return _amp / _d2**(2.5)
			else:
				raise ValueError("Zero or negative denominator in Plummer_Density")
		return Plummer_dens


# Plummer cumulative mass
def Plummer_Mass(mass = None, a = None):
	'''Returns the cumulative mass of a Plummer sphere at r'''
	if mass is None:
		raise ValueError("mass is a required parameter in Plummer_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in Plummer_Mass")
	elif a is None:
		raise ValueError("a is a required parameter in Plummer_Mass")
	elif a < 0:
		raise ValueError("a must be non-negative in Plummer_Mass")
	else:
		def Plummer_M(*r):
			_amp = mass
			_r2 = norm2(*r)
			_r = norm(*r)
			if _r > 0:
				return _amp * _r**3 / (_r2+a**2)**1.5
			else:
				raise ValueError("Zero or negative radius in Plummer_Mass")
		return Plummer_M


# Plummer velocity dispersion
def Plummer_VelDisp(mass = None, a = None):
	'''Returns the 1D Plummer radial velocity dispersion at r, sigma_r(r)^2
	See https://en.wikipedia.org/wiki/Plummer_model
	'''
	if mass is None:
		raise ValueError("mass is a required argument in Plummer_VelDisp")
	elif mass <= 0:
		raise ValueError("mass must be positive in Plummer_VelDisp")
	elif a is None:
		raise ValueError("a is a required argument in Plummer_VelDisp")
	elif a < 0:
		raise ValueError("a must be non-negative in Plummer_VelDisp")
	else:
		def Plummer_veldisp(*r):
			_amp = Grav*mass
			_r = norm(*r)
			_d = math.sqrt(_r**2+a**2)
			if _d > 0:
				return _amp / (6. * _d)
			else:
				raise ValueError("Zero or negative denominator in Plummer_VelDisp")
		return Plummer_veldisp



# NAVARRO, FRENK, AND WHITE (1997) MODEL

# NFW potential
def NFW_Potential(mass = None, rs = None):
	'''Implements a wrapper for the NFW (1997) potential'''
	if mass is None:
		raise ValueError("mass is a required argument in NFW_Potential")
	elif mass <= 0:
		raise ValueError("mass must be positive in NFW_Potential")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Potential")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Potential")
	else:

		def NFW_Pot(*r):
# 			_amp = -4. * math.pi * Grav * rho0 * rs**3
			_amp = -1.*Grav * mass
			_r = norm(*r)
			if _r > 0:
				return math.log(1.+_r/rs)*_amp/_r
			else:
				raise ValueError("Zero or negative radius in NFW_Potential")

		# private attributes (= parent func. params) to allow for access from outside
		NFW_Pot._rs = rs
		NFW_Pot._mass = mass
		return NFW_Pot


# NFW density profile
def NFW_Density(mass = None, rs = None):
	'''Returns the NFW density at r'''
	if mass is None:
		raise ValueError("mass is a required argument in NFW_Density")
	elif mass <= 0:
		raise ValueError("mass must be positive in NFW_Density")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Density")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Density")
	else:
		def NFW_dens(*r):
			rho0 = mass / (4.*math.pi*rs**3)
			_r = norm(*r)
			_x = _r / rs
			if _x > 0:
				return rho0 / (_x * (1.+_x)**2)
			else:
				raise ValueError("Zero or negative radius in NFW_Density")
		return NFW_dens


# NFW cumulative mass
def NFW_Mass(mass = None, rs = None):
	'''Returns the cumulative NFW mass at r'''
	if mass is None:
		raise ValueError("mass is a required argument in NFW_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in NFW_Mass")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Mass")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Mass")
	else:
		def NFW_M(*r):
# 			_amp = 4. * math.pi * rho0 * rs**3
			_amp = mass
			_r = norm(*r)
			return _amp * (math.log(1.+_r/rs) - _r/(_r+rs))
		return NFW_M


# NFW circular velocity
def NFW_Vcirc(mass = None, rs = None):
	'''Returns the spherical NFW circular velocity at r'''
	if mass is None:
		raise ValueError("mass is a required argument in NFW_Vcirc")
	elif mass <= 0:
		raise ValueError("mass must be positive in NFW_Vcirc")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Vcirc")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Vcirc")
	else:
		def NFW_Vc(*r):
			_r = norm(*r)
			if _r > 0:
				return math.sqrt((1./_r)*Grav*NFW_Mass(mass,rs)(*r))
			else:
				raise ValueError("Zero or negative radius in NFW_Vc")
		return NFW_Vc


# NFW maximum circular velocity
def NFW_Vmax(mass = None, rs = None):
	'''Returns the (approximate) maximum spherical NFW circular
	velocity roughly reached at 2.16*rs.'''
	if mass is None:
		raise ValueError("mass is a required argument in NFW_Vmax")
	elif mass <= 0:
		raise ValueError("mass must be positive in NFW_Vmax")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Vmax")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Vmax")
	else:
		return NFW_Vcirc(mass,rs)(2.16258*rs)


# NFW velocity dispersion
def NFW_VelDisp(mass = None, rs = None):
	'''Returns the approximate 1D NFW radial velocity dispersion at r, sigma_r(r)^2. The approximation is taken from Zentner and Bullock (2003,
	their equation 6).
	'''
	if mass is None:
		raise ValueError("mass is a required argument in NFW_VelDisp")
	elif mass <= 0:
		raise ValueError("mass must be positive in NFW_VelDisp")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_VelDisp")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_VelDisp")
	else:
		def NFW_veldisp(*r):
			_r = norm(*r)
			_x = _r / rs
			_vmax = NFW_Vmax(mass,rs)
			return (_vmax * (1.4393*_x**0.354/(1.+1.1756*_x**0.725)))**2
		return NFW_veldisp



# HERNQUIST (1990) MODEL

# Hernquist potential
def Hernquist_Potential(mass = None, a = None):
	'''Implements a wrapper for the Hernquist (1990) potential.
	Reduces to Kepler potential for a=0.'''
	if mass is None:
		raise ValueError("mass is a required parameter in Hernquist_Potential")
	elif mass <= 0:
		raise ValueError("mass must be positive in Hernquist_Potential")
	elif a is None:
		raise ValueError("a is a required parameter in Hernquist_Potential")
	elif a < 0:
		raise ValueError("a must be non-negative in Hernquist_Potential")
	else:
		def Hernquist_Pot(*r):
			_amp = Grav*mass
			_r = norm(*r)
			_d = _r + a
			if _d > 0:
				return -1. * (_amp / _d)
			else:
				raise ValueError("Zero or negative denominator in Hernquist_Potential")
		# private attributes (= parent func. params) to allow for access from outside
		Hernquist_Pot._mass = mass
		Hernquist_Pot._a = a
		return Hernquist_Pot


# Hernquist density profile
def Hernquist_Density(mass = None, a = None):
	'''Returns the Hernquist density at r'''
	if mass is None:
		raise ValueError("mass is a required parameter in Hernquist_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in Hernquist_Mass")
	elif a is None:
		raise ValueError("a is a required parameter in Hernquist_Mass")
	elif a < 0:
		raise ValueError("a must be non-negative in Hernquist_Mass")
	else:
		def Hernquist_dens(*r):
			_amp = mass * a / (2.*math.pi)
			_r = norm(*r)
			_d = _r + a
			if _r > 0 and _d > 0:
				return _amp / (_r * _d**3)
			else:
				raise ValueError("Zero or negative denominator in Hernquist_Density")
		return Hernquist_dens


# Hernquist cumulative mass
def Hernquist_Mass(mass = None, a = None):
	'''Returns the cumulative Hernquist mass at r'''
	if mass is None:
		raise ValueError("mass is a required parameter in Hernquist_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in Hernquist_Mass")
	elif a is None:
		raise ValueError("a is a required parameter in Hernquist_Mass")
	elif a < 0:
		raise ValueError("a must be non-negative in Hernquist_Mass")
	else:
		def Hernquist_M(*r):
			_r2 = norm2(*r)
			_r = norm(*r)
			_d = _r + a
			if _d > 0:
				return mass * _r2 / _d**2
			else:
				raise ValueError("Zero or negative denominator in Hernquist_Mass")
		return Hernquist_M


# Hernquist velocity dispersion
def Hernquist_VelDisp(mass = None, a = None):
	'''Returns the 1D Hernquist radial velocity dispersion at r, sigma_r(r)^2,
	as given by Hernquist (1990, his equation 10)
	See http://adsabs.harvard.edu/abs/1990ApJ...356..359H
	'''
	if mass is None:
		raise ValueError("mass is a required parameter in Hernquist_VelDisp")
	elif mass <= 0:
		raise ValueError("mass must be positive in Hernquist_VelDisp")
	elif a is None:
		raise ValueError("a is a required parameter in Hernquist_VelDisp")
	elif a <= 0:
		raise ValueError("a must be positive in Hernquist_VelDisp")
	else:
		def Hernquist_veldisp(*r):
			_amp = Grav*mass/(12.*a)
			_r = norm(*r)
			if _r > 0:
				_rpa = _r + a
				_roa = _r / a
				_term1 = 12.*_r*_rpa**3*math.log(_rpa/_r)/a**4
				_term2 = _r/_rpa
				_term3 = 25. + _roa*(52 + _roa*(42. + _roa*12.))
				if _term1 >= _term2*_term3:
					return _amp*(_term1 - _term2*_term3)
				else:
					raise ValueError("Negative velocity dispersion in Hernquist_VelDisp")
			else:
				raise ValueError("Zero or negative radius in Hernquist_VelDisp")
		return Hernquist_veldisp



# PSEUDO-ISOTHERMAL SPHERE (PITS)
# Used in e.g.:
# Jiang & Binney (2000, https://ui.adsabs.harvard.edu/abs/2000MNRAS.314..468J/abstract) -> with exp. tapper
# Kormendy & Freeman (2016, https://ui.adsabs.harvard.edu/abs/2016ApJ...817...84K/abstract)

# PITS potential
def PITS_Potential(mass = None, a = None):
	'''Implements a wrapper for the PITS potential'''
	if mass is None:
		raise ValueError("mass is a required argument in PITS_Potential")
	elif mass <= 0:
		raise ValueError("mass must be positive in PITS_Potential")
	elif a is None:
		raise ValueError("a is a required argument in PITS_Potential")
	elif a <= 0:
		raise ValueError("a must be positive in PITS_Potential")
	else:

		def PITS_Pot(*r):
# 			_amp = 4. * math.pi * Grav * rho0 * a**2
			_amp = Grav * mass / a
			_r = norm(*r)
			_x = _r / a
			if _x > 0:
				return _amp * ( 0.5 * math.log(1.+_x**2) + math.atan(_x)/_x )
			else:
				raise ValueError("Zero or negative radius in PITS_Potential")

		# private attributes (= parent func. params) to allow for access from outside
		PITS_Pot._a = a
		PITS_Pot._mass = mass
		return PITS_Pot


# PITS density profile
def PITS_Density(mass = None, a = None):
	'''Returns the PITS density at r'''
	if mass is None:
		raise ValueError("mass is a required argument in PITS_Density")
	elif mass <= 0:
		raise ValueError("mass must be positive in PITS_Density")
	elif a is None:
		raise ValueError("a is a required argument in PITS_Density")
	elif a <= 0:
		raise ValueError("a must be positive in PITS_Density")
	else:
		def PITS_dens(*r):
			rho0 = mass / (4.*math.pi*a**3)
			_r2 = norm2(*r)
			_x2 = _r2 / a**2
			return rho0 / (1. + _x2)
		return PITS_dens


# PITS cumulative mass
def PITS_Mass(mass = None, a = None):
	'''Returns the cumulative PITS mass at r'''
	if mass is None:
		raise ValueError("mass is a required argument in PITS_Mass")
	elif mass <= 0:
		raise ValueError("mass must be positive in PITS_Mass")
	elif a is None:
		raise ValueError("a is a required argument in PITS_Mass")
	elif a <= 0:
		raise ValueError("a must be positive in PITS_Mass")
	else:
		def PITS_M(*r):
# 			_amp = 4. * math.pi * rho0 * a**3
			_amp = mass
			_r = norm(*r)
			_x = _r / a
			return _amp * (_x - math.atan(_x))
		return PITS_M


# PITS asymptotic circular velocity (V_infinity)
# See de Blok (2010), their equation 1 and below:
# https://ui.adsabs.harvard.edu/abs/2010AdAst2010E...5D/abstract)
def PITS_Vinf(mass = None, a = None):
	'''Returns the asymptotic PITS circular velocity'''
	if mass is None:
		raise ValueError("mass is a required argument in PITS_Vinf")
	elif mass <= 0:
		raise ValueError("mass must be positive in PITS_Vinf")
	elif a is None:
		raise ValueError("a is a required argument in PITS_Vinf")
	elif a <= 0:
		raise ValueError("a must be positive in PITS_Vinf")
	else:
# 		return math.sqrt(4. * math.pi * Grav * rho0 * a**2)
		return math.sqrt(Grav * mass / a)


# PITS velocity dispersion
def PITS_VelDisp(mass = None, a = None):
	'''Returns the 1D PITS  radial velocity dispersion at r, sigma_r(r)^2.
	The result is taken from Kormendy & Freeman (2016, their
	equation 4). Beware that they may have adopted a different
	convention for v_inf as we do here (see Kormendy & Freeman (2016,
	https://ui.adsabs.harvard.edu/abs/2016ApJ...817...84K/abstract;
	their equation 1 and bottom of page.
	'''
	if mass is None:
		raise ValueError("mass is a required argument in PITS_VelDisp")
	elif mass <= 0:
		raise ValueError("mass must be positive in PITS_VelDisp")
	elif a is None:
		raise ValueError("a is a required argument in PITS_VelDisp")
	elif a <= 0:
		raise ValueError("a must be positive in PITS_VelDisp")
	else:
		def PITS_veldisp(*r):
			_r = norm(*r)
			_x = _r / a
			if _x > 0:
				_atanx = math.atan(_x)
				_vinf2 = (PITS_Vinf(mass,a))**2
				return _vinf2 * (1.+_x**2) * (0.125*math.pi**2 - _atanx / _x - 0.5 * _atanx**2)
			else:
				raise ValueError("Zero or negative radius in PITS_VelDisp")
		return PITS_veldisp




# FIELDS AND FORCES

def tidal_radius(*r, m1_func = None, m2_func = None):
	'''Returns the tidal radius of a body b2 -- described by its mass
	distribution m2_func -- orbiting at a distance r from a central mass b1,
	described by its mass distribution m1_func.
	It requires a root finding algorithm as the tidal radius is implicitly
	defined through a relation between the cumulative masses of the bodies
	which depend on the bodys' relative position and the tidal radius.
	See Klypin et al. (1999a, their equation 6) or
	Nichols & Bland-Hawthorn (2009, their equation 12)
	IMPORTANT: The following assumptions are made:
	- both masses are spherically symmetric
	- the tidal radius is not expected to be smaller than 1.e-5*r
	'''
	if m1_func is None:
		raise ValueError("m1_func is a required parameter in tidal_radius")
	elif m2_func is None:
		raise ValueError("m2_func is a required parameter in tidal_radius")
	elif len(r) < 1:
		raise ValueError("r must have dimension >= 1 in tidal_radius, but has {}".format(len(r)))
	else:
		_r = norm(*r)
		dMdr = grad_r(*r, func = m1_func)
		m1 = m1_func(_r)
		def func(_rt):
			m2 = m2_func(_rt)
			return _rt**3*(2.-(_r/m1)*dMdr)*m1 - _r**3*m2
		# find a range where root is bracketed
		r0, r1 = 1.e-5*_r, _r
		while func(r0)*func(r1)>0 and r0<r1:
			r0 *= 1.1
		return brent_root(f = func, x0 = r0, x1 = r1, max_iter=50, tolerance=1.e-5)


def tidal_radius_approx(*r, m1_func = None, m2_func = None):
	'''Returns the approximate tidal radius of a body b2 -- described by its mass
	distribution m2_func -- orbiting at a distance r from a central mass b1,
	described by its mass distribution m1_func.
	It requires a root finding algorithm as the tidal radius is implicitly
	defined through a relation between the cumulative masses of the bodies
	which depend on the bodys' relative position and the tidal radius.
	See Jiang & Loeb (2000, their equation 5)
	See Dierickx & Loeb (2017a, their equation 8)
	Note that their equations differ by a numerical factor of order 1.
	IMPORTANT: The following assumptions are made:
	- both masses are spherically symmetric
	- the tidal radius is not expected to be smaller than 1.e-5*r
	'''
	if m1_func is None:
		raise ValueError("m1_func is a required parameter in tidal_radius_approx")
	elif m2_func is None:
		raise ValueError("m2_func is a required parameter in tidal_radius_approx")
	elif len(r) < 1:
		raise ValueError("r must have dimension >= 1 in tidal_radius_approx, but has {}".format(len(r)))
	else:
		_r = norm(*r)
		def func(_rt):
			_rrel = max(1.e-5*_r,abs(_r - _rt))
			return 2.*_rt**3 * m1_func(_rrel) - _rrel**3 * m2_func(_rt)
		# find a range where root is bracketed
		r0, r1 = 1.e-5*_r, _r
		while func(r0)*func(r1) > 0:
			r0 *= 1.1
		return brent_root(f = func, x0 = 1.e-5*_r, x1 = _r, max_iter=50, tolerance=1.e-5)


def mass_bound(m1_func = None, m2_func = None):
	'''Returns the mass of an object b2 -- described by its mass
	distribution m2_func -- within its tidal radius orbiting at a distance
	r from another object b1, described by its mass distribution m1_func.
	IMPORTANT: Spherically symmetry of the masses is assumed!
	'''
	if m1_func is None:
		raise ValueError("m1_func is a required parameter in mass_bound")
	elif m2_func is None:
		raise ValueError("m2_func is a required parameter in mass_bound")
	else:
		# it must be a generic function of t and r
		def _mass_b(t=None,r=None):
			_rt = [tidal_radius(*r,m1_func=m1_func,m2_func=m2_func)]
			mb = m2_func(*_rt)
			return mb
		return _mass_b



def dyn_friction_simpl():
	'''Calcualates the *magnitude* of the deceleration experienced by an object of
	mass M due to dynamical friction exerted by a surrounding, uniform density
	field of matter composed of particles of mass m using a simplified expression
	for the sole purpose of testing.
	r - position vector
	v - velocity vector
	mass - mass M
	rho - density field (function of r)
	veldisp - unused parameter; introduced for consistency with other dyn.frict.funcs.
	See https://en.wikipedia.org/wiki/Dynamical_friction
	'''
	def dynfric_simpl(r = None, v = None, mass = None, rho = None, veldisp = None):
		if r is None:
			raise ValueError("r is a required parameter in dyn_friction_simpl")
		elif v is None:
			raise ValueError("v is a required parameter in dyn_friction_simpl")
		elif mass is None:
			raise ValueError("mass is a required parameter in dyn_friction_simpl")
		elif rho is None:
			raise ValueError("rho is a required parameter in dyn_friction_simpl")
		elif len(r) != len(v):
			raise ValueError("r and v must have the same dimensionality in dyn_friction_simpl")
		else:
			const = -4. * math.pi # arbitrary
			amp = const * Grav**2 * mass
			_v = norm(*v)
			if _v > 0:
				return amp * rho(*r) / _v**3	# v^3 because it is the magnitude
	return dynfric_simpl


def dyn_friction_maxwell(eps = None):
	'''Calcualates the negative acceleration (or deceleration) of an object of
	mass M due to dynamical friction exerted by a surrounding, uniform density
	field of matter (host) composed of particles of mass m using the Chandrasekhar
	formula. Assumptions:
	1 - M >> m
	2 - The velocity of the matter particles obeys a Maxwellian distribution.
	See https://en.wikipedia.org/wiki/Dynamical_friction
	Notes:
	- eps is the 'softening length' of the satellite and is a tunable
	parameter (see e.g. Hashimoto et al. 2003)
	- it is not yet clear which velocity dispersion is required (1D or 3D)
	'''
	if eps is None:
		raise ValueError("eps is a required parameter in dyn_friction_maxwell")
	else:

		def dynfric_maxwell(r = None, v = None, mass = None, rho = None, veldisp = None):
			if r is None:
				raise ValueError("r is a required parameter in dyn_friction_maxwell")
			elif v is None:
				raise ValueError("v is a required parameter in dyn_friction_maxwell")
			elif mass is None:
				raise ValueError("mass is a required parameter in dyn_friction_maxwell")
			elif rho is None:
				raise ValueError("rho is a required parameter in dyn_friction_maxwell")
			elif veldisp is None:
				raise ValueError("veldisp is a required parameter in dyn_friction_maxwell")
			elif len(r) != len(v):
				raise ValueError("r and v must have the same dimensionality in dyn_friction_maxwell")
			else:
				const = -4. * math.pi
				_log_lambda = couloumb_log_hfm03(*r,epsilon = eps)
				_v = norm(*v)
				if _v > 0:
					amp = const * Grav**2 * mass * _log_lambda * rho(*r) / _v**3
# 					sigma = math.sqrt(vel_disp_1d_jeans(*r))	# too expensive!
					sigma = math.sqrt(veldisp(*r))
					if sigma > 0:
						_X = _v / (math.sqrt(2.) * sigma)
						erfX = math.erf(_X)
						XexpX2 = 2.*_X*math.exp(-1.*_X**2)/math.sqrt(math.pi)
						return amp * (erfX - XexpX2)
					else:
						raise ValueError("Vanishing velocity dispersion in dyn_friction_maxwell")
				else:
					raise ValueError("v must be positive in dyn_friction_maxwell")
		return dynfric_maxwell


def vel_disp_1d_jeans(*r, rho = None, pot = None):
	'''
	WARNING: This routine is computationally expensive.
	It is better to use the analytic expression specific for
	each potential.

	returns the 1-dimensional velocity dispersion calculated
	from the spherically symmetric Jeans (1915) equation,
	valid for a non-rotating, isotropic (beta=0) system.
	The equation is integrated from r to 'infinity'.
	The gradient of the potential is calculated using a central
	finite difference scheme of order 4 and constant step size.
	The integration is performed using a Runge-Kutta method of order
	4 with constant step size
	'''
	if rho is None:
		raise ValueError("rho is a required parameter in vel_disp_1d_jeans")
	elif pot is None:
		raise ValueError("pot is a required parameter in vel_disp_1d_jeans")
	elif len(r) < 1:
		raise \
		ValueError("r must be at least of dimension 1 in vel_disp_1d_jeans, but is {}".format(len(r)))
	else:
		_r = norm(*r)
		if _r >= Infinity:
			raise ValueError("Upper integration limit ('infinity') exceeded in vel_disp_1d_jeans")
		else:
			def integrand(x,*y):
				dpotdr = grad_r(x,func=pot)
				return rho(x) * dpotdr
			# note: the parameter 'order' refers to order of the ODE, *not* of the method
			_steps = 1000
			_h = (Infinity-_r) / _steps
			ics = [_r, 0.]
			_, int = ode_rk4(dydx=[integrand], order=1, initVal=ics, steps=_steps, stepSize=_h)
			return int[0][_steps]/rho(*r)


def couloumb_log_hfm03(*r, epsilon = None):
	'''Returns the Couloumb logarithm used to calculate the dynamical friction
	deceleration, using the parametrization by Hashimoto et al. (2003, their
	equation 5), appropriate for N-body systems.
	'''
	if epsilon is None:
		raise ValueError("epsilon is a required parameter in couloumb_log_hfm03")
	elif epsilon <= 0:
		raise ValueError("epsilon must be positive in couloumb_log_hfm03")
	else:
		_r = norm(*r)
		denom = 1.4 * epsilon
		if _r <= 0:
			raise ValueError("r must be positive in couloumb_log_hfm03")
		elif _r < denom:	# avoid dynamical friction 'acceleration'
			return 0.
		else:
			return math.log(_r / denom)


def grav_field(*r, pot = None):
	'''Returns the gravitational field g = -Grad Phi at r for a given
	potential Phi. Grad Phi is calculated using a central finite difference
	scheme of order 4 and constant step of size 1.e-4.'''
	if pot is None:
		raise ValueError("pot is a required parameter in grav_field")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got {} in grav_field".format(len(r)))
	else:
		return [ 1.*f for f in grad(*r,func=pot)]


def v_circ(*r, amp = None):
	'''Returns the circular velocity at r for a Kepler potential'''
	if amp is None:
		raise ValueError("amp is a required parameter in v_circ")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got {} in v_circ".format(len(r)))
	else:
		_r = norm(*r)
		if _r > 0:
			return math.sqrt(amp / _r)
		else:
			raise ValueError("Zero or negative r in v_circ")


def v_circ_gen(*r, pot = None):
	'''Returns the circular velocity at r for any potential;
	given by sqrt(r Grad Phi).
	Grad Phi is calculated using a central finite difference
	scheme of order 4 and constant step of size 1.e-4.'''
	if pot is None:
		raise ValueError("pot is a required parameter in v_circ_gen")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got {} in v_circ_gen".format(len(r)))
	else:
		r_dot_grad_pot = 0.
		for i in range(len(r)):
			r_dot_grad_pot += r[i]*cen_diff_first(*r, var=i, func=pot, delta_x=1.e-4, order=4)
		return math.sqrt(r_dot_grad_pot)



# ENERGIES

def ePot(*r,pot = None):
	'''Specific potential energy'''
	if pot is None:
		raise ValueError("pot is a required parameter in ePot")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got {} in ePot".format(len(r)))
	else:
		return pot(*r)


def eKin(*v):
	'''Specific kinetic energy'''
	if len(v) < 1:
		raise ValueError("expected at least one element in v, got {} in eKin".format(len(v)))
	else:
		return 0.5 * norm2(*v)


# ORBITAL FUNCTIONS

def vis_viva(*r, amp = None, a = None):
	'''Returns the speed from the Vis-viva equation.
	Only valid for elliptic orbits.'''
	_r = norm(*r)
	if 2*a > _r:
		return math.sqrt(amp*(2./_r - 1./a))
	else:
		raise ValueError("Undefined velocity in vis-viva; set a > r/2")


# Semimajor axis from vis-viva
def semimajor_from_vis(r = None,v = None, amp = None):
	'''Returns from vis-viva equation; only valid for elliptic orbits.
	Superseded by function semimajor().'''
	if r is None:
		raise ValueError("r is a required parameter in semimajor_from_vis")
	elif v is None:
		raise ValueError("v is a required parameter in semimajor_from_vis")
	elif amp is None:
		raise ValueError("amp is a required parameter in semimajor_from_vis")
	else:
		v2 = v**2
		k1 = 2./r
		k2 = v2/amp
		if k1 != k2:
			return 1./(k1 - k2)
		else:
			print("seminajor axis undefined in semimajor_from_vis")
			return float('NaN')
