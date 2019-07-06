'''
Author:	Thorsten Tepper Garcia
Date:	26/06/2019
'''

import math
import sys
sys.path.insert(0,"../.")
from config.phys_consts import Grav
from num_diff.central_diff import cen_diff_first


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
		raise ValueError("r is a required parameter in cross_prod.")
	elif v is None:
		raise ValueError("v is a required parameter in cross_prod.")
	elif len(r) != len(v):
		raise ValueError("r and v must have the same length in cross_prod.")
	elif len(r) < 2:
		raise ValueError("cannot calculate the cross product in less than 2D in cross_prod.")
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
		raise ValueError("r is a required parameter in dot_prod.")
	elif v is None:
		raise ValueError("v is a required parameter in dot_prod.")
	elif len(r) != len(v):
		raise ValueError("r and v must have the same length in dot_prod.")
	elif len(r) < 2:
		raise ValueError("cannot calculate the dot product in less than 2D in dot_prod.")
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
		raise ValueError("vec is a required parameter in rodrigues_rot.")
	elif axis is None:
		raise ValueError("axis is a required parameter in rodrigues_rot.")
	elif incl is None:
		raise ValueError("incl is a required parameter in rodrigues_rot.")
	elif len(vec) != len(axis):
		raise ValueError("vec and axis must have the same length in rodrigues_rot.")
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



# POTENTIALS AND CORRESPONDING MASS PROFILES

def Kepler_Potential(amp = None):
	'''Implements a wrapper for the Kepler potential where amp = G*M.'''
	if amp is None:
		raise ValueError("amp is a required argument in Kepler_Potential.")
	elif amp <= 0:
		raise ValueError("amp must be positive in Kepler_Potential.")
	else:
		def Kepler_Pot(*r):
			_amp = amp
			_r = norm(*r)
			if _r > 0:
				return -1. * (_amp / _r)
			else:
				raise ValueError("Zero or negative radius in Kepler_Potential.")
		return Kepler_Pot

# Kepler mass
# Note that this function is trivial, but is included here
# and used in the code for completeness and self-consistency.
def Kepler_Mass(mass = None):
	if mass is None:
		raise ValueError("mass is a required parameter in Kepler_Mass.")
	else:
		def Kepler_M(*r):
			return mass
		return Kepler_M


# Plummer potential
def Plummer_Potential(amp = None, a = None):
	'''Implements  a wrapper for the Plummer (1911) potential
	where amp = G*M. Reduces to Kepler potential for a=0.'''
	if amp is None:
		raise ValueError("amp is a required argument in Plummer_Potential.")
	elif amp <= 0:
		raise ValueError("amp must be positive in Plummer_Potential.")
	elif a is None:
		raise ValueError("a is a required argument in Plummer_Potential.")
	elif a < 0:
		raise ValueError("a must be non-negative in Plummer_Potential.")
	else:
		def Plummer_Pot(*r):
			_r2 = norm2(*r)
			_r = math.sqrt(_r2+a**2)
			if _r > 0:
				return -1. * (amp / _r)
			else:
				raise ValueError("Zero or negative radius in Plummer_Potential.")
		# private attributes (= parent func. params) to allow for access from outside
		Plummer_Pot._a = a
		return Plummer_Pot


# Plummer cumulative mass
def Plummer_Mass(mass = None, a = None):
	if mass is None:
		raise ValueError("mass is a required parameter in Plummer_Mass")
	elif a is None:
		raise ValueError("a is a required parameter in Plummer_Mass")
	elif a < 0:
		raise ValueError("a must be non-negative in Plummer_Mass")
	else:
		def Plummer_M(*r):
			amp = mass
			_r2 = norm2(*r)
			_r = norm(*r)
			if _r > 0:
				return amp * _r**3 / (_r2+a**2)**1.5
			else:
				raise ValueError("Zero or negative radius in Plummer_Mass.")
		return Plummer_M


# Navarro, Frenk, and White (1997) potential
def NFW_Potential(rho0 = None, rs = None):
	'''Implements  a wrapper for the NFW (1997) potential'''
	if rho0 is None:
		raise ValueError("rho0 is a required argument in NFW_Potential.")
	elif rho0 <= 0:
		raise ValueError("rho0 must be positive in NFW_Potential.")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Potential.")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Potential.")
	else:

		def NFW_Pot(*r):
			amp = -4. * math.pi * Grav * rho0 * rs**3
			_r = norm(*r)
			if _r > 0:
				return math.log(1.+_r/rs)*amp/_r
			else:
				raise ValueError("Zero or negative radius in NFW_Potential.")

		# private attributes (= parent func. params) to allow for access from outside
		NFW_Pot._rs = rs
		NFW_Pot._rho0 = rho0
		return NFW_Pot


# NFW cumulative mass
def NFW_Mass(rho0 = None, rs = None):
	'''Returns the cumulative NFW mass at r'''
	if rho0 is None:
		raise ValueError("rho0 is a required argument in NFW_Mass.")
	elif rho0 <= 0:
		raise ValueError("rho0 must be positive in NFW_Mass.")
	elif rs is None:
		raise ValueError("rs is a required argument in NFW_Mass.")
	elif rs <= 0:
		raise ValueError("rs must be positive in NFW_Mass.")
	else:
		def NFW_M(*r):
			amp = 4. * math.pi * rho0 * rs**3
			_r = norm(*r)
			return amp * (math.log(1.+_r/rs) - _r/(_r+rs))
		return NFW_M


# Hernquist potential
def Hernquist_Potential(amp = None, a = None):
	'''Implements a wrapper for the Hernquist (1990) potential
	where amp = G*M. Reduces to Kepler potential for a=0.'''
	if amp is None:
		raise ValueError("amp is a required parameter in Hernquist_Potential")
	elif a is None:
		raise ValueError("a is a required parameter in Hernquist_Potential")
	elif a < 0:
		raise ValueError("a must be non-negative in Hernquist_Potential")
	else:
		def Hernquist_Pot(*r):
			_r = norm(*r)
			if _r > 0:
				return -1. * (amp / (_r+a))
			else:
				raise ValueError("Zero or negative radius in Hernquist_Potential.")
		# private attributes (= parent func. params) to allow for access from outside
		Hernquist_Pot._a = a
		return Hernquist_Pot
	

# Hernquist cumulative mass
def Hernquist_Mass(mass = None, a = None):
	if mass is None:
		raise ValueError("mass is a required parameter in Hernquist_Mass")
	elif a is None:
		raise ValueError("a is a required parameter in Hernquist_Mass")
	elif a < 0:
		raise ValueError("a must be non-negative in Hernquist_Mass")
	else:
		def Hernquist_M(*r):
			_r2 = norm2(*r)
			_r = norm(*r)
			if _r > 0:
				return mass * _r2 / (_r+a)**2
			else:
				raise ValueError("Zero or negative radius in Hernquist_Mass.")
		return Hernquist_M


# FIELDS AND FORCES

def grav_field(*r, pot = None):
	'''Returns the gravitational field g = Grad Phi at r for a given
	potential Phi. Grad Phi is calculated using a central finite difference
	scheme of order 4.'''
	if pot is None:
		raise ValueError("pot is a required parameter in grav_field.")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got{} in grav_field.".format(len(r)))
	else:
		dpotdr = []
		for i in range(len(r)):
			dpotdr.append(cen_diff_first(*r, var=i, func=pot, delta_x=1.e-4, order=4))
		return dpotdr


def v_circ(*r, amp = None):
	'''Returns the circular velocity at r for a Kepler potential'''
	if amp is None:
		raise ValueError("amp is a required parameter in v_circ.")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got{} in v_circ.".format(len(r)))
	else:
		_r = norm(*r)
		if _r > 0:
			return math.sqrt(amp / _r)
		else:
			raise ValueError("Zero or negative r in v_circ.")

def v_circ_gen(*r, pot = None):
	'''Returns the circular velocity at r for any potential;
	given by sqrt(r Grad Phi).
	Grad Phi is calculated using a central finite difference
	scheme of order 4.'''
	if pot is None:
		raise ValueError("pot is a required parameter in v_circ_gen.")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got{} in v_circ_gen.".format(len(r)))
	else:
		r_dot_dpotdr = 0.
		for i in range(len(r)):
			r_dot_dpotdr += r[i]*cen_diff_first(*r, var=i, func=pot, delta_x=1.e-4, order=4)
		return math.sqrt(r_dot_dpotdr)



# ENERGIES

def ePot(*r,pot = None):
	'''Specific potential energy'''
	if pot is None:
		raise ValueError("pot is a required parameter in ePot.")
	elif len(r) < 1:
		raise ValueError("expected at least one element in r, got{} in ePot.".format(len(r)))
	else:
		return pot(*r)


def eKin(*v):
	'''Specific kinetic energy'''
	if len(v) < 1:
		raise ValueError("expected at least one element in v, got{} in eKin.".format(len(v)))
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
		raise ValueError("Undefined velocity in vis-viva; set a > r/2.")


# Semimajor axis from vis-viva
def semimajor_from_vis(r = None,v = None, amp = None):
	'''Returns from vis-viva equation; only valid for elliptic orbits.
	Superseded by function semimajor().'''
	if r is None:
		raise ValueError("r is a required parameter in semimajor_from_vis.")
	elif v is None:
		raise ValueError("v is a required parameter in semimajor_from_vis.")
	elif amp is None:
		raise ValueError("amp is a required parameter in semimajor_from_vis.")
	else:
		v2 = v**2
		k1 = 2./r
		k2 = v2/amp
		if k1 != k2:
			return 1./(k1 - k2)
		else:
			print("seminajor axis undefined in semimajor_from_vis.")
			return float('NaN')