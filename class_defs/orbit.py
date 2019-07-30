'''
Author:	Thorsten Tepper Garcia
'''

import math
import sys
from utils import funcs
import config.phys_consts as pc
from num_diff.central_diff import cen_diff_first	# central finite difference scheme (recommended)
													# alternative: forward finite difference scheme
from ode_int.leapfrog import ode_leap				# leapfrog time integrator (recommended)
													# alternative: euler-richardson
from config import units


class Orbit():
	"""General Orbit class representing an orbit"""
	def __init__(self, b1 = None, b2 = None):
		"""
		NAME:

			__init__

		PURPOSE:

			Initialize an Orbit instance

		INPUT:

			b1 - a Body instance
			b2 - a Body instance

		OUTPUT:

			Instance          

		HISTORY:

			2019-07-05 - Written - TTG

		"""
		if b1 is None or b2 is None:
			raise ValueError("Two Body instances are required when initialising an Orbit object.")
		else:
			self.b1 = b1
			self.b2 = b2
			self.orbit_info()	# dump initial orbital parameters to stdout


	def pos(self):
		'''position of body 2 relative to body 1'''
		return self.b2.pos_rel(self.b1)


	def vel(self):
		'''velocity of body 2 relative to body 1'''
		return self.b2.vel_rel(self.b1)


	def dist(self):
		'''distance of body 2 relative to body 1'''
		return self.b2.dist_rel(self.b1)


	def speed(self):
		'''speed of body 2 relative to body 1'''
		return self.b2.speed_rel(self.b1)


# Deactivated because not consistent
# 	def mtot(self):
# 		'''total mass'''
# 		return self.b1.mass_scale + self.b2.mass_scale


# Deactivated because not consistent
# 	def mred(self):
# 		'''reduced mass'''
# 		return self.b1.mass_scale*self.b2.mass_scale/self.mtot()


	def grav_param(self):
		'''gravitational parameter'''
		_r = self.pos()
		_mtot_r = self.b1.mass_cum(*_r)+self.b2.mass_cum(*_r)
		return pc.Grav*_mtot_r


	def ang_mom_vec(self):
		'''Angular momentum of body 2 relative to body 1'''
		return funcs.cross_prod(self.pos(),self.vel())


	def ang_mom(self):
		'''Magnitude ofngular momentum of body 2 relative to body 1'''
		return funcs.norm(*self.ang_mom_vec())


	def ang_mom_vec_norm(self):
		'''Normalised angular momentum of body 2 relative to body 1'''
		_L = self.ang_mom()
		if _L > 0:
			return [ l/_L for l in self.ang_mom_vec()]
		else:
			raise ValueError("Initial conditions imply a purely radial orbit (vanishing angular momentum).")


	def v_rad(self):
		'''radial velocity'''
		_r = self.pos()
		_d = self.dist()
		_v = self.vel()
		if _d > 0:
			return funcs.dot_prod(_r,_v) / _d
		else:
			raise ValueError("Non-positive distance in v_rad()")


	def v_tan(self):
		'''tangential velocity'''
		_r = self.pos()
		_d = self.dist()
		_v = self.vel()
		if _d > 0:
			return funcs.norm(*funcs.cross_prod(_r,_v)) / _d
		else:
			raise ValueError("Non-positive distance in v_tan()")


	def energy_pot(self):
		'''specific potential energy'''
		_r = self.pos()
		return (self.b1.potential(*_r)+self.b2.potential(*_r))


	def energy_kin(self):
		'''specific relative kinetic energy'''
		_v = self.speed()
		return 0.5*_v**2


	def energy_tot(self):
		'''total specific energy'''
		return self.energy_pot()+self.energy_kin()


	def eccentricity_vec(self):
		'''calculates the specific Laplace-Runge-Lenz vector a.k.a. 'eccentricity vector'
		'specific' means it is normalised by (G Mtot)*(M_reduced)
		the direction of the LRL vector lies along the symmetry axis of the
		conic section and points from the center of force toward the periapsis,
		i.e., the point of closest approach.
		See: https://en.wikipedia.org/wiki/Laplace–Runge–Lenz_vector'''
		_vxh = funcs.cross_prod(self.vel(),self.ang_mom_vec())
		_r = self.pos()
		_d = self.dist()
		_mu = self.grav_param()
		return [a/_mu - b/_d for a, b in zip(_vxh, _r)]


	def eccentricity(self):
		'''scalar eccentricity'''
		return funcs.norm(*self.eccentricity_vec())


	def semi_latus_rec(self):
		'''semi-latus rectum'''
		return (self.ang_mom())**2/self.grav_param()


	def semimajor(self):
		'''semimajor axis; undefined for e = 1'''
		e = self.eccentricity()
		if e < 1 or e > 1:
			return self.semi_latus_rec() / (1. - e**2)
		else:
			print("\nWARNING:\n\tsemimajor axis not defined for e = 1.\n")
			return float('NaN')


	def semiminor(self):
		'''semiminor axis; undefined for e = 1'''
		e = self.eccentricity()
		if e < 1 or e > 1:
			return self.semi_latus_rec() / math.sqrt(abs(1. - e**2))
		else:
			print("\nWARNING:\n\tsemiminor axis not defined for e = 1.\n")
			return float('NaN')


	def apocenter(self):
		'''apocentric distance; undefined for e >= 1'''
		e = self.eccentricity()
		if e < 1:
			return	self.semi_latus_rec() / (1. - e)
		else:
			print("\nWARNING:\n\tapocentre not defined for e >= 1.\n")
			return float('NaN')


	def pericenter(self):
		'''pericentric distance'''
		return	self.semi_latus_rec() / (1. + self.eccentricity())


	def v_apo(self):
		'''speed at apocenter'''
		return self.ang_mom() / self.apocenter()


	def v_peri(self):
		'''speed at pericenter'''
		return self.ang_mom() / self.pericenter()


	def circumference(self):
		''' Returns the approximate cirumference of an orbit with 0 < e < 1,
		exact for e = 0, and undefined for e >= 1.'''
		a = self.semimajor()
		b = self.semiminor()
		if a>=b:
			# correction factor
			corr = 0.25 * ((a-b)/(a+b))**2
			# approximate circumference to order 4 accuracy;
			# higher orders expansion exist!
			return math.pi * (a+b) * (1.+corr+(corr/4.)**4)
		else:
			print("\nWARNING:\n\torbital circumference not defined for e >= 1.\n")
			return float('NaN')


	def period(self):
		'''orbital period; undefined for e >= 1'''
		a = self.semimajor()
		if a > 0:
			return 2.*math.pi*a*math.sqrt(a/self.grav_param())
		else:
			print("\nWARNING:\n\torbital period not defined for a < 0.\n")
			return float('NaN')


	def pericentric_freq(self):
		'''orbital frequency at pericenter; undefined for e >= 1'''
		return self.v_peri() / self.circumference()


	def energy_drift_lim(self):
		'''returns an estimate of the time step required to
		avoid energy drift in the time integration.'''
		w = self.pericentric_freq()
		if not math.isnan(w):
			return 1. / w / (math.sqrt(2.)*math.pi)
		else:
			return float('NaN')


	def orbital_plane(self, vec = None):
		'''Determine transformation that maps state vectors onto
		orbital plane; makes use of Rodrigues' rotation formula
		See: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
		The transformation is such that the transformed angular
		momentum aligns with the positive z-axis.
		If vec is not provided, the function returns the inclination
		angle (radian) and the reference unit vector.'''
		_zaxis = [0.,0.,1.]
		_l = self.ang_mom_vec_norm()
		_cos_theta = funcs.dot_prod(_zaxis,_l)
		incl = math.acos(_cos_theta)
		# order of cross producto very important!
		_k = funcs.cross_prod(_l,_zaxis)
		k_norm = funcs.norm(*_k)
		if k_norm > 0:
			k_vec = [k/k_norm for k in _k]
		else:
			k_vec = _k
		if vec is None:
			return incl, k_vec
		else:
			return funcs.rodrigues_rot(vec,k_vec,incl)


	def apsidal_angle(self):
		'''Calculates the angle between the periapsis -- given by
		the eccentricity vector -- and the x-axis on the orbital plane.'''
		_evec = self.eccentricity_vec()
		_evecrot = self.orbital_plane(_evec)
		_enorm = funcs.norm(*_evecrot)
		if _enorm > 1.e-10:
			return math.atan2(_evecrot[1],_evecrot[0])
		else:
			return 0.


	def asc_node_vec(self):
		'''vector of the ascending node;
		the reference vector is taken to be the z-axis;
		the reference plane, the xy-plane.'''
		_k = [0.,0.,1.]
		_h = self.ang_mom_vec()
		return funcs.cross_prod(_k,_h)


	def long_asc_node(self):
		'''Calculate the longitude of the ascending node in radians.
		The reference vector is taken to be the z-axis; the reference
		plane, the xy-plane.
		See https://en.wikipedia.org/wiki/Longitude_of_the_ascending_node'''
		_vec = self.asc_node_vec()
		_norm = funcs.norm(*_vec)
		if _norm > 0:
			_n = [n/_norm for n in _vec]
			if _n[1] >= 0:
				return math.acos(_n[0])
			else:
				return 2.*math.pi - math.acos(_n[0])
		else:
			print("\nWARNING:\n\tlongitude of ascending node undefined. Set to 0 by convention.\n")
			return 0.


	def arg_periapsis(self):
		'''Calculates the argument of periapsis in radians.
		The reference vector is taken to be the z-axis;
		the reference plane, the xy-plane.'''
		_n = self.asc_node_vec()
		_nnorm = funcs.norm(*_n)
		_e = self.eccentricity_vec()
		_enorm = self.eccentricity()
		if _nnorm*_enorm > 0:
			_ndote = funcs.dot_prod(_n,_e)
			if _e[2] >= 0:
				return math.acos(_ndote/_nnorm/_enorm)
			else:
				return 2.*math.pi - math.acos(_ndote/_nnorm/_enorm)
		elif _enorm > 1.e-10:
			_l = self.ang_mom_vec()
			if _l[2] >= 0:
				return math.atan2(_e[1],_e[0])
			else:
				return 2.*math.pi - math.atan2(_e[1],_e[0])
		else:
			return 0.


	def kepler_orbit_polar(self,theta):
		'''Evaluates the analytic kepler orbit:
			r(theta) = p / 1 + e cos(theta-theta0)
			at theta (given in radians).'''
		_p = self.semi_latus_rec()
		_e = self.eccentricity()
		_theta0 = self.apsidal_angle()
		return _p / (1.+ _e * math.cos(theta-_theta0))


	def kepler_orbit_cartesian(self,theta):
		'''Returns the cartesian x,y coordinates of a Kepler orbit on the orbital
		plane'''
		_r = self.kepler_orbit_polar(theta)
		return _r*math.cos(theta), _r*math.sin(theta)



	def orbit_info(self, filename = None):
		'''prints a number of relevant orbit parameter values to stdout or a file if provided'''

		_incl, _ = self.orbital_plane()	# skip k_vec from output
		_orbital_period_peri = 1./self.pericentric_freq()

		if filename is not None:
			fileout = open(filename,'wt')
			first_chr = '# '
		else:
			fileout = sys.stdout
			first_chr = ''

		print("\n{}Initial (osculating) orbital parameters of the system:\n".format(first_chr), \
			file = fileout)
		print("{}{:>40}{:>15}".format(first_chr,"Potential of body 1:",self.b1.potential.__name__), \
			file = fileout)
		print("{}{:>40}{:>15}".format(first_chr,"Potential of body 2:",self.b2.potential.__name__), \
			file = fileout)
		print("{}{:>40}{:15.4E}".format(first_chr,"Bound mass of body 1:",self.b1.mass_bound), \
			file = fileout)
		print("{}{:>40}{:15.4E}\n".format(first_chr,"Bound mass of body 2:",self.b2.mass_bound), \
			file = fileout)
# 		print("{}{:>40}{:15.4E}\n".format(first_chr,"Reduced mass:",self.mred()), \
# 			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Initial rel. distance (r21_0):",self.dist()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Initial rel. speed (v21_0):",self.speed()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Initial tangential vel. (v_tan_0):",self.v_tan()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Initial radial vel. (v_rad_0):",self.v_rad()), \
			file = fileout)
		print("{}{:>40} ({:5.3f},{:5.3f},{:5.3f})".format(first_chr,"Rel. specific ang. mom. vec. (h_vec):", \
			*self.ang_mom_vec()), file = fileout)
		print("{}{:>40} ({:5.3f},{:5.3f},{:5.3f})".format(first_chr,"[normalised]:", \
			*self.ang_mom_vec_norm()), file = fileout)
		print("{}{:>40} ({:5.3f},{:5.3f},{:5.3f})".format(first_chr,"[along z-axis]:", \
			*self.orbital_plane(self.ang_mom_vec())), file = fileout)
		print("{}{:>40}{:15.4f}\n".format(first_chr,"Rel. specific ang. mom. (h):",self.ang_mom()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Eccentricity (e):",self.eccentricity()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Semi-latus rectum (p):",self.semi_latus_rec()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Semimajor axis (a):",self.semimajor()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Semiminor axis (b):",self.semiminor()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Orbital inclination (psi_0; deg):",math.degrees(_incl)), \
			file = fileout)
		print("{}{:>40} ({:5.3f},{:5.3f},{:5.3f})".format(first_chr,"Ascending node (n_vec):", \
			*self.asc_node_vec()), file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Long. of asc. node (Omega_0; deg):", \
			math.degrees(self.long_asc_node())), file = fileout)
		print("{}{:>40} ({:5.3f},{:5.3f},{:5.3f})".format(first_chr,"Eccentricity vector (e_vec):", \
			*self.eccentricity_vec()), file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Apsidal angle (phi_0; deg):", \
			math.degrees(self.apsidal_angle())), file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Argument of periapsis (omega_0; deg):", \
			math.degrees(self.arg_periapsis())), file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Rel. pericentre (rp):",self.pericenter()), \
			file = fileout)
		print("{}{:>40}{:15.4f}".format(first_chr,"Vel. at pericentre (vp):",self.v_peri()), \
			file = fileout)
		if self.eccentricity() < 1.:
			print("{}{:>40}{:15.4f}".format(first_chr,"Rel. apocentre (ra):",self.apocenter()), \
			file = fileout)
			print("{}{:>40}{:15.4f}".format(first_chr,"Vel. at apocentre (va):",self.v_apo()), \
			file = fileout)
			print("{}{:>40}{:15.4f}".format(first_chr,"Rel. orbital period (T):",self.period()), \
			file = fileout)
			print("{}{:>40}{:15.4f}".format(first_chr,"Approx. orbit circumference (u):",self.circumference()), \
			file = fileout)
			print("{}{:>40}{:15.4f}".format(first_chr,"Approx. pericentric period (Tp):",_orbital_period_peri), \
			file = fileout)
		print("{}{:>40}{:15.4E}".format(first_chr,"Rel. specific potential energy (V):",self.energy_pot()), \
			file = fileout)
		print("{}{:>40}{:15.4E}".format(first_chr,"Rel. specific kinetic energy (T):",self.energy_kin()), \
			file = fileout)
		print("{}{:>40}{:15.4E}\n".format(first_chr,"Rel. specific total energy (E):",self.energy_tot()), \
			file = fileout)
		
		if filename is not None:
			fileout.close()



	def integrate(self, time_start = None, time_end = None, time_step = None):
		"""
		NAME:

			integrate

		PURPOSE:

			Integrate an orbit in time

		INPUT:

			time_start - initial time
			time_end - end time
			time_step - size of time step

		OUTPUT:
          
          time_steps - Number of time steps
          time - list of integration time steps (with length = time_steps)
          EoM - interleaved list if state vectors for each body (x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
        		for each time step

		HISTORY:

			2019-07-09 - Written - TTG
			2019-07-11 - Moved declaration of accelerations to independent function - TTG

		"""
		if time_start is None:
			raise ValueError("time_start is a required parameter in integrate.")
		elif time_end is None:
			raise ValueError("time_end is a required parameter in integrate.")
		elif time_step is None:
			raise ValueError("time_step is a required parameter in integrate.")
		else:

			# Set up time integrator
			backwards_switch = 1
			if time_end < 0:
				backwards_switch = -1
				time_end = backwards_switch * time_end
				
			# number of time steps
			_T = (time_end - time_start)
			_N = math.ceil( _T / time_step)

			# initial conditions (flip velocity sign if integrating backwards)
			# note: the shape of _ics is dictated by module ode_leap
			_ics = \
				[time_start, \
				self.b1.x, backwards_switch * self.b1.vx, \
				self.b1.y, backwards_switch * self.b1.vy, \
				self.b1.z, backwards_switch * self.b1.vz, \
				self.b2.x, backwards_switch * self.b2.vx, \
				self.b2.y, backwards_switch * self.b2.vy, \
				self.b2.z, backwards_switch * self.b2.vz]

			# accelerations:
			# dynamical friction is turned into a repulsive force when backwards_switch < 0
			# note: the shape of _F is dictated by time integration module (ode_leap or ode_eulric)
			_F = [*self.set_accelerations(backwards_switch)]

			print("\nOrbit integration{}".format(" (backwards):" if backwards_switch < 0 else " (forward):"))
			
			if backwards_switch > 0:
				print("\tTime range [t0,t1] = [{},{}]".format(time_start,time_end))
			else:
				print("\tTime range [t0,t1] = [{},{}]".format(backwards_switch*time_end,time_start))
			print("\tTime steps {};  Step size: {:8.2E}".format(_N,time_step))
			
			if self.b1.dynamical_friction is None and self.b2.dynamical_friction is None:
				print("\tMaximum time step to avoid energy drift: {:12.6E}".format(self.energy_drift_lim()))
				# see http://physics.bu.edu/py502/lectures3/cmotion.pdf:
				print("\tExpected accumulated error in x and v of order {:8.2E}".format((_T * time_step)**2))

			# Integrate system
			# Note: x1 = EoM[0], vx1 = EoM[1], y1 = EoM[2], vy1 = EoM[3], z1 = EoM[4], vz1 = EoM[5],
			#       x2 = EoM[6], vx2 = EoM[7], y2 = EoM[8], vy2 = EoM[9], z2 = EoM[10], vz2 = EoM[11]
			# 		each EoM[i] is 'function' of time, i.e. an array where each item corresponds to a
			# 		given integration time, i.e. EoM[i][t]
			time, EoM = ode_leap(dr2dt2 = _F, rank = 12, initCond = _ics, steps = _N, stepSize = time_step)

			
			# invert time arrow and velocities after backwards integration
			if backwards_switch < 1:
				for i in range(len(time)):
					time[i] = backwards_switch * time[i]
					for j in range(1,len(EoM),2):
						EoM[j][i] = backwards_switch * EoM[j][i]
				time.reverse()
				for j in range(len(EoM)):
						EoM[j].reverse()

			return time, EoM


	def set_accelerations(self, bwi_switch=None):
		"""
		NAME:

			set_accelerations

		PURPOSE:

			Set acceleration function in equation of motion for orbit integration.
			The functions defined below correspond each to one of the (numerical)
			partial derivatives of the potential of each body.
			Here, a central finite difference scheme is used to calculate the derivatives.
			Other schemes are possible, but a CFD has been tested thoroughly and has been found
			to perform well, even for cuspy potentials such as the Kepler potential.

		INPUT:

			bwi_switch - backwards integration switch; turns dissipative forces into repulsive
			forces for backwards integration

		OUTPUT:

			dvx1dt - acceleration function of body 1 along x
			dvy1dt - acceleration function of body 1 along y
			dvz1dt - acceleration function of body 1 along z
			dvx2dt - acceleration function of body 2 along x
			dvy2dt - acceleration function of body 2 along y
			dvz2dt - acceleration function of body 2 along z

		HISTORY:

			2019-07-11 - Written - TTG
		
		NOTES:
		
			Could be integrated instead into the Body class (TTG)

		"""
		if bwi_switch is None:
			bwi_switch = 1
		
		# Central finite difference scheme parameters:
		# Notes:
		# - a 2nd order scheme conserves well both energy and angular momentum.
		# - a 4th order scheme conserves energy well, and angular momentum generally to machine precision.
		# DO NOT change the values of the following parameters unless absolutely necessary!
		intStep = 1.e-4
		accOrder = 4

		# Note: the array f = (x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
		# the shape of f is dictated by the module cen_diff_first
		# Recall: Force field = - Grad Phi
		
		# Dynamical friction: can only be exerted by one body onto the other!
		if self.b1.dynamical_friction is not None:

			def dvx1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=0, func=self.b2.potential, delta_x=intStep, order=accOrder)
				return force_grav

			def dvy1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=1, func=self.b2.potential, delta_x=intStep, order=accOrder)
				return force_grav

			def dvz1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=2, func=self.b2.potential, delta_x=intStep, order=accOrder)
				return force_grav

			def dvx2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				_vx21 = f[7] - f[1]
				_vy21 = f[9] - f[3]
				_vz21 = f[11] - f[5]
				_vvec = [_vx21,_vy21,_vz21]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=0, func=self.b1.potential, delta_x=intStep, order=accOrder)
				force_df = bwi_switch * \
					self.b1.dynamical_friction(r=_rvec, v=_vvec, mass=self.b2.mass_evol_x(t,_rvec,_vvec,bwi_switch), rho=self.b1.dens, \
						veldisp=self.b1.vel_disp) * _vx21
				return force_grav+force_df

			def dvy2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				_vx21 = f[7] - f[1]
				_vy21 = f[9] - f[3]
				_vz21 = f[11] - f[5]
				_vvec = [_vx21,_vy21,_vz21]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=1, func=self.b1.potential, delta_x=intStep, order=accOrder)
				force_df = bwi_switch * \
					self.b1.dynamical_friction(r=_rvec, v=_vvec, mass=self.b2.mass_evol_y(t,_rvec,_vvec,bwi_switch), rho=self.b1.dens, \
						veldisp=self.b1.vel_disp) * _vy21
				return force_grav+force_df

			def dvz2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				_vx21 = f[7] - f[1]
				_vy21 = f[9] - f[3]
				_vz21 = f[11] - f[5]
				_vvec = [_vx21,_vy21,_vz21]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=2, func=self.b1.potential, delta_x=intStep, order=accOrder)
				force_df = bwi_switch * \
					self.b1.dynamical_friction(r=_rvec, v=_vvec, mass=self.b2.mass_evol_z(t,_rvec,_vvec,bwi_switch), rho=self.b1.dens, \
						veldisp=self.b1.vel_disp) * _vz21
				return force_grav+force_df
		
		elif self.b2.dynamical_friction is not None:

			def dvx1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				_vx12 = f[1]-f[7]
				_vy12 = f[3]-f[9]
				_vz12 = f[5]-f[11]
				_vvec = [_vx12,_vy12,_vz12]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=0, func=self.b2.potential, delta_x=intStep, order=accOrder)
				force_df = bwi_switch * \
					self.b2.dynamical_friction(r=_rvec, v=_vvec, mass=self.b1.mass_evol_x(t,_rvec,_vvec,bwi_switch), rho=self.b2.dens, \
						veldisp=self.b2.vel_disp) * _vx12
				return force_grav+force_df

			def dvy1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				_vx12 = f[1]-f[7]
				_vy12 = f[3]-f[9]
				_vz12 = f[5]-f[11]
				_vvec = [_vx12,_vy12,_vz12]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=1, func=self.b2.potential, delta_x=intStep, order=accOrder)
				force_df = bwi_switch * \
					self.b2.dynamical_friction(r=_rvec, v=_vvec, mass=self.b1.mass_evol_y(t,_rvec,_vvec,bwi_switch), rho=self.b2.dens, \
						veldisp=self.b2.vel_disp) * _vy12
				return force_grav+force_df

			def dvz1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				_vx12 = f[1]-f[7]
				_vy12 = f[3]-f[9]
				_vz12 = f[5]-f[11]
				_vvec = [_vx12,_vy12,_vz12]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=2, func=self.b2.potential, delta_x=intStep, order=accOrder)
				force_df = bwi_switch * \
					self.b2.dynamical_friction(r=_rvec, v=_vvec, mass=self.b1.mass_evol_z(t,_rvec,_vvec,bwi_switch), rho=self.b2.dens, \
						veldisp=self.b2.vel_disp) * _vz12
				return force_grav+force_df

			def dvx2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=0, func=self.b1.potential, delta_x=intStep, order=accOrder)
				return force_grav

			def dvy2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=1, func=self.b1.potential, delta_x=intStep, order=accOrder)
				return force_grav

			def dvz2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				force_grav = \
					-1.*cen_diff_first(*_rvec, var=2, func=self.b1.potential, delta_x=intStep, order=accOrder)
				return force_grav
		
		else:	# no dynamical friction

			def dvx1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				return  -1.*cen_diff_first(*_rvec, var=0, func=self.b2.potential, delta_x=intStep, order=accOrder)

			def dvy1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				return   -1.*cen_diff_first(*_rvec, var=1, func=self.b2.potential, delta_x=intStep, order=accOrder)

			def dvz1dt(t,*f):
				_x12 = f[0]-f[6]
				_y12 = f[2]-f[8]
				_z12 = f[4]-f[10]
				_rvec = [_x12,_y12,_z12]
				return   -1.*cen_diff_first(*_rvec, var=2, func=self.b2.potential, delta_x=intStep, order=accOrder)

			def dvx2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				return  -1.*cen_diff_first(*_rvec, var=0, func=self.b1.potential, delta_x=intStep, order=accOrder)

			def dvy2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				return  -1.*cen_diff_first(*_rvec, var=1, func=self.b1.potential, delta_x=intStep, order=accOrder)

			def dvz2dt(t,*f):
				_x21 = f[6]-f[0]
				_y21 = f[8]-f[2]
				_z21 = f[10]-f[4]
				_rvec = [_x21,_y21,_z21]
				return  -1.*cen_diff_first(*_rvec, var=2, func=self.b1.potential, delta_x=intStep, order=accOrder)

		return dvx1dt, dvy1dt, dvz1dt, dvx2dt, dvy2dt, dvz2dt


	def write_table(self, time_list = None, state_vector = None, filename = None, output_freq = None):
		"""
		NAME:

			write_table

		PURPOSE:

			Dump the orbital evolution to an ascii table.
			The table consists of a total of 26 columns. The first 16 contain, in that
			order, the time, the specific relative angular momentum, the relative potential
			energy, the relative kinetic energy, and (x,vx,y,vy,z,vz) for each of the
			bodies. The next 8 columnns contain (x',vx',y',vy') for each of the
			bodies, where the primed coordinates and velocities correspond to the coordinates
			and velocities on the orbital plane. In other words, these coordinates represent a
			rotated version of the intrinsic orbit such that the relative angular momentum
			of the system aligns with the z-axis. If the intrinsic orbit has this property,
			the primed and unprimed coordinates are identical. Note that the primed z
			coordinates are ignored, because they all vanish by definition. The last two
			columns contain the (x,y) coordinates of the analytic solution of the problem,
			using the orbital parameters self-consistently calculated within the code.
			The latter can be directly compared to the set of primed coordinates as a check
			of the numerical result against the expected analytic solution, keeping in
			mind that for non-Keplerian potentials the analytic solution merely corresponds
			to an osculating orbit consistent with the initial conditions.


		INPUT:

			filename - full path to output file
			time_list - list of integration time steps (with length = time_steps)
			state_vector - interleaved list if state vectors for each body:
				  					(x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
				  				for each time step
		OUTPUT:
          
			None

		HISTORY:

			2019-07-09 - Written - TTG
			2019-07-23 - Move cons.laws to independent routine - TTG

		"""
		if time_list is None:
			raise ValueError("time_list is a required parameter in write_table")
		elif state_vector is None:
			raise ValueError("state_vector is a required parameter in write_table")
		elif filename is None:
			raise ValueError("filename is a required parameter in write_table")
		elif output_freq is None:
			raise ValueError("output_freq is a required parameter in write_table")
		else:

			# include initial orbital parameters in header of output file
			self.orbit_info(filename)

			# append orbital evolution to file
			f = open(filename, 'at')

			f.write(("{:<9}{:<6}{:8} {:<14}"+"{:6} {:<9}"*2+"{:4} {:<9}"*22+"\n").\
				format("# time", units.TIME_UNIT_STR, \
					"ang.mom.", units.ANG_MOM_UNIT_STR, \
					"ePot", units.ENERGY_UNIT_STR, "eKin", units.ENERGY_UNIT_STR, \
					"x1", units.LENGTH_UNIT_STR, "vx1", units.VEL_UNIT_STR, \
					"y1", units.LENGTH_UNIT_STR, "vy1", units.VEL_UNIT_STR, \
					"z1", units.LENGTH_UNIT_STR, "vz1", units.VEL_UNIT_STR, \
					"x2", units.LENGTH_UNIT_STR, "vx2", units.VEL_UNIT_STR, \
					"y2", units.LENGTH_UNIT_STR, "vy2", units.VEL_UNIT_STR, \
					"z2", units.LENGTH_UNIT_STR, "vz2", units.VEL_UNIT_STR, \
					"x1_proj", units.LENGTH_UNIT_STR, "vx1_proj", units.VEL_UNIT_STR, \
					"y1_proj", units.LENGTH_UNIT_STR, "vy1_proj", units.VEL_UNIT_STR, \
					"x2_proj", units.LENGTH_UNIT_STR, "vx2_proj", units.VEL_UNIT_STR, \
					"y2_proj", units.LENGTH_UNIT_STR, "vy2_proj", units.VEL_UNIT_STR, \
					"KeplerOrbitAna_X", units.LENGTH_UNIT_STR, \
					"KeplerOrbitAna_Y", units.LENGTH_UNIT_STR))

			time_steps = len(time_list)

			for t in range(0,time_steps,output_freq):

				x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 = \
					state_vector[0][t], state_vector[1][t], state_vector[2][t], state_vector[3][t], \
					state_vector[4][t], state_vector[5][t], state_vector[6][t], state_vector[7][t], \
					state_vector[8][t], state_vector[9][t], state_vector[10][t], state_vector[11][t]
				r1 = [x1,y1,z1]
				v1 = [vx1,vy1,vz1]
				r2 = [x2,y2,z2]
				v2 = [vx2,vy2,vz2]
				r_rel = [x2-x1,y2-y1,z2-z1]
				v_rel = [vx2-vx1,vy2-vy1,vz2-vz1]
				Ltot_vec = funcs.cross_prod(r_rel,v_rel)		# rel. spec. ang. mom. (divided by Mred)
				Ltot = funcs.norm(*Ltot_vec)
				ePot = \
					(self.b1.potential(*r_rel) + \
					self.b2.potential(*r_rel))						# rel. spec. potential energy (?)
				eKin = funcs.eKin(*v_rel)							# rel. spec. kin. energy
				r1_rot = self.orbital_plane(r1)					# rotate vectors onto orbital plane
				v1_rot = self.orbital_plane(v1)
				r2_rot = self.orbital_plane(r2)
				v2_rot = self.orbital_plane(v2)
				x1_rot,y1_rot,_ = r1_rot							# ignore z-comps. because they are 0.
				vx1_rot,vy1_rot,_ = v1_rot
				x2_rot,y2_rot,_ = r2_rot
				vx2_rot,vy2_rot,_ = v2_rot

				f.write(("{:<15.8f}{:<23.10E}"+"{:<16.8E}"*2+"{:<14.4E}"*22+"\n").\
					format(time_list[t]*units.TIME.value, Ltot, ePot, eKin, \
								x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2, \
								x1_rot, vx1_rot, y1_rot, vy1_rot, x2_rot, vx2_rot, y2_rot, vy2_rot, \
								*self.kepler_orbit_cartesian(t*2.*math.pi/time_steps)\
							)\
						 )

			f.close()

			print("\nOutput with timestep frequency {} written to file:\n\n\t{}\n".\
				format(output_freq,filename))

			# conservation laws
			self.conservation(time=time_list,sv=state_vector)
			
			# Improve on this output; e.g. make it dependent on whether mass
			# evolution was explicitly set by input parameter.
			# mass evolution (if set)
			self.ini_end_out(time=time_list,sv=state_vector)

			return 0


	def ini_end_out(self, time = None, sv = None):
		"""
		NAME:

			ini_end_out

		PURPOSE:

			Dump to stodut some useful information of the initial
			and final state vectors and masses. The latter is only
			relevant if some form of mass evolution has been switched
			on


		INPUT:

			time - list of integration time steps (with length = time_steps)
			sv - interleaved list if state vectors for each body:
				  	(x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
				  for each time step

		OUTPUT:
          
			None

		HISTORY:

			2019-07-23 - Written - TTG

		"""
		if time is None:
			raise ValueError("time is a required parameter in write_table")
		elif sv is None:
			raise ValueError("sv is a required parameter in ini_end_out")
		else:
			# set index for initial and final time steps
			t_ini = 0
			t_end = len(time)-1
			print("\nInitial / final state vectors:")
			x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 = \
					sv[0][t_ini], sv[1][t_ini], sv[2][t_ini], sv[3][t_ini], \
					sv[4][t_ini], sv[5][t_ini], sv[6][t_ini], sv[7][t_ini], \
					sv[8][t_ini], sv[9][t_ini], sv[10][t_ini], sv[11][t_ini]
			r_ini = [x2-x1,y2-y1,z2-z1]
			v_ini = [vx2-vx1,vy2-vy1,vz2-vz1]
			x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 = \
					sv[0][t_end], sv[1][t_end], sv[2][t_end], sv[3][t_end], \
					sv[4][t_end], sv[5][t_end], sv[6][t_end], sv[7][t_end], \
					sv[8][t_end], sv[9][t_end], sv[10][t_end], sv[11][t_end]
			r_end = [x2-x1,y2-y1,z2-z1]
			v_end = [vx2-vx1,vy2-vy1,vz2-vz1]
			print("\tInfall relative position: ({:.5f},{:.5f},{:.5f})".format(*r_ini))
			print("\tInfall relative distance: {:.5f}".format(funcs.norm(*r_ini)))
			print("\tInfall relative velocity: ({:.5f},{:.5f},{:.5f})".format(*v_ini))
			print("\tInfall relative speed: {:.5f}".format(funcs.norm(*v_ini)))
			print("\tPresent-day relative position: ({:.5f},{:.5f},{:.5f})".format(*r_end))
			print("\tPresent-day relative distance: {:.5f}".format(funcs.norm(*r_end)))
			print("\tPresent-day relative velocity: ({:.5f},{:.5f},{:.5f})".format(*v_end))
			print("\tPresent-day relative speed: {:.5f}\n".format(funcs.norm(*v_end)))
			# This may not yet be entirely consistent:
			print("The following are only relevant is masses are allowed to evolved:")
			if time[0] < 0:
				out_string = "Infall"
			else:
				out_string = "Present-day"
			print("\t{} bound mass of body 1: {:E}".format(out_string,self.b1.mass_bound))
			print("\t{} bound mass of body 2: {:E}\n".format(out_string,self.b2.mass_bound))

			return 0



	def conservation(self, time = None, sv = None):
		"""
		NAME:

			conservation

		PURPOSE:

			Calculate the state of conservation laws, i.e. the change in the
			total relative specific energy, and the angular momentum vector (both
			its direction and magnitude)


		INPUT:

			time - list of integration time steps (with length = time_steps)
			sv - interleaved list if state vectors for each body:
				 	(x1,vx1,y1,vy1,z1,vz1,x2,vx2,y2,vy2,z2,vz2)
				  for each time step

		OUTPUT:
          
			None

		HISTORY:

			2019-07-23 - Written - TTG

		"""
		if time is None:
			raise ValueError("time is a required parameter in write_table")
		elif sv is None:
			raise ValueError("sv is a required parameter in write_table")
		else:

			# Initial (reference) values
			ePot0 = self.energy_pot()
			eKin0 = self.energy_kin()
			Ltot_vec0 = self.ang_mom_vec()
			Ltot0 = self.ang_mom()
			x_axis = [1.,0.,0.]						# used to check direction of L
			y_axis = [0.,1.,0.]						# used to check direction of L
			z_axis = [0.,0.,1.]						# used to check direction of L

			eTot_Cons, lCons, lvecCons = 0., 0., 0.
			cos_ax0 = funcs.dot_prod(x_axis,Ltot_vec0) / Ltot0
			cos_ay0 = funcs.dot_prod(y_axis,Ltot_vec0) / Ltot0
			cos_az0 = funcs.dot_prod(z_axis,Ltot_vec0) / Ltot0
			time_steps = len(time)

			for t in range(0,time_steps):

				x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2 = \
					sv[0][t], sv[1][t], sv[2][t], sv[3][t], \
					sv[4][t], sv[5][t], sv[6][t], sv[7][t], \
					sv[8][t], sv[9][t], sv[10][t], sv[11][t]
				r1 = [x1,y1,z1]
				v1 = [vx1,vy1,vz1]
				r2 = [x2,y2,z2]
				v2 = [vx2,vy2,vz2]
				r_rel = [x2-x1,y2-y1,z2-z1]
				v_rel = [vx2-vx1,vy2-vy1,vz2-vz1]
				Ltot_vec = funcs.cross_prod(r_rel,v_rel)		# rel. spec. ang. mom. (divided by Mred)
				Ltot = funcs.norm(*Ltot_vec)
				ePot = \
					(self.b1.potential(*r_rel) + \
					self.b2.potential(*r_rel))						# rel. spec. potential energy (?)
				eKin = funcs.eKin(*v_rel)							# rel. spec. kin. energy

				# Calculate state of conservation laws
				eTot_t = abs((ePot+eKin)/(ePot0+eKin0)-1.)
				if eTot_t > eTot_Cons:
					eTot_Cons = eTot_t
				lCons_t = abs((Ltot/Ltot0)-1.)
				if lCons_t > lCons:
					lCons = lCons_t
				cos_ax = funcs.dot_prod(x_axis,Ltot_vec) / Ltot
				cos_ay = funcs.dot_prod(y_axis,Ltot_vec) / Ltot
				cos_az = funcs.dot_prod(z_axis,Ltot_vec) / Ltot
				if abs(cos_ax0) > 0.:
					lvecConsx_t = abs(cos_ax/cos_ax0-1.)
				else:
					lvecConsx_t = 0.
				if abs(cos_ay0) > 0.:
					lvecConsy_t = abs(cos_ay/cos_ay0-1.)
				else:
					lvecConsy_t = 0.
				if abs(cos_az0) > 0.:
					lvecConsz_t = abs(cos_az/cos_az0-1.)
				else:
					lvecConsz_t = 0.
				lvecCons_t = max(lvecConsx_t,max(lvecConsy_t,lvecConsz_t))
				if lvecCons_t > lvecCons:
					lvecCons = lvecCons_t 

			print("Conservation laws:")
			if self.b1.dynamical_friction is not None or self.b2.dynamical_friction is not None:
				print("\n\nWARNING:\n\tDynamical friction is switched on!")
				print("\tNeither energy nor angular momentum will be conserved.\n")
			print("{:>60} {:8.3E} %.".\
				format("Energy conservation to better than",1.e2*eTot_Cons))
			print("{:>60} {:8.3E} %.".\
				format("Angular momentum conservation (magnitude) to better than",1.e2*lCons))
			print("{:>60} {:8.3E} %.\n".\
				format("Angular momentum conservation (direction) to better than",1.e2*lvecCons))

			if eTot_t > 1. or lCons_t > 1:
				print("\nWARNING:\n\tPossible merger scenario! Inspect orbit carefully.\n")

			return 0

