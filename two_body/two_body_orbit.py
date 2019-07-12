#!/opt/local/bin//python3.5

'''
Author:			Thorsten Tepper Garcia
Last modified:	09/07/2019

See README for information on the code's background, usage, etc.

'''

# Import necessary modules

# built-in modules
import sys

# costum modules
sys.path.insert(0,"../.")							# include top directory in module search path
from utils import io
from class_defs import orbit


# Gather initial conditions
ics = io.get_input()

# Initialise Orbit object
orbit = orbit.Orbit(ics.body1,ics.body2)

# Integrate orbit
time, EoM = orbit.integrate(ics.t_0,ics.t_1,ics.delta_t)

# Write orbit evolution to file
outDir = "./output/"
outFile = outDir + ics.filename_prefix + "_out.dat"

orbit.write_table(time_list = time, state_vector = EoM, filename = outFile, output_freq = 10)
