'''
Author:	Thorsten Tepper Garcia
Date:	02/07/2019
'''

import sys
sys.path.insert(0,"./init")							# include initial conditions directory
import importlib									# needed to import ICs' as module
from class_defs import body

def get_input():
	"""
	NAME:

		get_input

	PURPOSE:

		Gather initial conditions for and initialize an body instances

	INPUT:

		None

	OUTPUT:

		ics - Object with initial conditons parameters

	HISTORY:

		2019-07-11 - Written - TTG

	"""
	if len(sys.argv) < 2:
		print("\nUSAGE:")
		print("{} <input parameter file>\n".format(sys.argv[0]))
		exit()
	else:
		relative_path='../two_body/init'
		ics_file = sys.argv[1]
		ext = ".py"
		if ext in ics_file:
			print("\nRemoved unnecessary extension {} in input parameter file.".format(ext))
			ics_file = ics_file.replace(ext,'') 	# remove file extension if present

		# Set initial conditions
		print("\nGathering input parameters from file {}...".format(ics_file))
		ic = importlib.import_module(ics_file,package=relative_path)

		# Time integration settings
		try:
			ic.t_0
		except:
			ic.t_0 = 0.
		try:
			ic.t_1
		except:
			ic.t_1 = 10.
		try:
			ic.delta_t
		except:
			ic.delta_t = 0.001

		# sanity check in case of backwards integration
		if ic.t_1 < 0 and ic.t_0 != 0:
			raise \
			ValueError("Initial time must be 0 for backwards integration")

		# Body 1
		try:
			ic.Mass1
		except:
			raise ValueError("Body mass is a required input parameter")
		try:
			r10_vec = [ic.x1_0,ic.y1_0,ic.z1_0]
		except:
			raise ValueError("Body coordinates are a required input parameter")
		try:
			v10_vec = [ic.vx1_0,ic.vy1_0,ic.vz1_0]
		except:
			raise ValueError("Body velocity is a required input parameter")
		try:
			Pot1 = ic.Potential1
		except:
			raise ValueError("Body potential is a required input parameter")
		try:
			DF1 = ic.Dynamical_Friction1
		except:
			DF1 = None
		try:
			mass1_evol = ics.Mass1_evol
		except:
			mass1_evol = None
		ic.body1 = body.Body(mass=ic.Mass1,pot=Pot1,r_vec=r10_vec,v_vec=v10_vec,df=DF1,massevol=mass1_evol)

		# Body 2
		try:
			ic.Mass2
		except:
			raise ValueError("Body mass is a required input parameter")
		try:
			r20_vec = [ic.x2_0,ic.y2_0,ic.z2_0]
		except:
			raise ValueError("Body coordinates are a required input parameter")
		try:
			v20_vec = [ic.vx2_0,ic.vy2_0,ic.vz2_0]
		except:
			raise ValueError("Body velocity is a required input parameter")
		try:
			Pot2 = ic.Potential2
		except:
			raise ValueError("Body potential is a required input parameter")
		try:
			DF2 = ic.Dynamical_Friction2
		except:
			DF2 = None
		try:
			mass2_evol = ic.Mass2_evol
		except:
			mass2_evol = None
		ic.body2 = body.Body(mass=ic.Mass2,pot=Pot2,r_vec=r20_vec,v_vec=v20_vec,df=DF2,massevol=mass2_evol)

		
		# sanity check
		if DF1 is not None and DF2 is not None:
			raise ValueError("Both bodies cannot excert dynamical friction simultaneously")
		else:
			print("Done.")

		# save init file name
		ic.filename_prefix = ics_file

		return ic

