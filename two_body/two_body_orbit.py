#!/opt/local/bin//python3.5

'''
Author:			Thorsten Tepper Garcia
Last modified:	09/07/2019

See README for information on the code's background, usage, etc.

'''

# Import necessary modules

# built-in modules
import importlib									# needed to import ICs' as module
import sys

# costum modules
sys.path.insert(0,"../.")							# include top directory in module search path
sys.path.insert(0,"./init")							# include initial conditions directory
from utils import io
from class_defs import body
from class_defs import orbit


# Collect input argument
initialConds = io.get_input()

# Set initial conditions
# (consider moving this into a function or class)
ic = importlib.import_module(initialConds)
print("\nGathering input parameters from file {}...".format(initialConds))

# Time integration settings
try:
	t0 = ic.t_0
except:
	t0 = 0.
try:
	t1 = ic.t_1
except:
	t1 = 10.
try:
	timeStep = ic.delta_t
except:
	timeStep = 0.001

# Body 1
try:
	Pot1 = ic.Potential1
except:
	Pot1 = None
try:
	DF1 = ic.Dynamical_Friction1
except:
	DF1 = None
r10_vec = [ic.x1_0,ic.y1_0,ic.z1_0]
v10_vec = [ic.vx1_0,ic.vy1_0,ic.vz1_0]
body1 = body.Body(mass=ic.Mass1,pot=Pot1,r_vec=r10_vec,v_vec=v10_vec,df=DF1)

# Body 2
try:
	Pot2 = ic.Potential2
except:
	Pot2 = None
try:
	DF2 = ic.Dynamical_Friction2
except:
	DF2 = None
r20_vec = [ic.x2_0,ic.y2_0,ic.z2_0]
v20_vec = [ic.vx2_0,ic.vy2_0,ic.vz2_0]
body2 = body.Body(mass=ic.Mass2,pot=Pot2,r_vec=r20_vec,v_vec=v20_vec,df=DF2)

# sanity check
if DF1 is not None and DF2 is not None:
	raise ValueError("Both bodies cannot excert dynamical friction simultaneously")
else:
	print("Done.")


# Initialise Orbit object
orbit = orbit.Orbit(body1,body2)

# Integrate orbit
time, EoM = orbit.integrate(t0,t1,timeStep)

# Write orbit evolution to file
outDir = "./output/"
outFile = outDir + initialConds + "_out.dat"

orbit.write_table(time_list = time, state_vector = EoM, filename = outFile, output_freq = 10)
