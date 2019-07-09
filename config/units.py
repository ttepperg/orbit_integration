'''
Author:	Thorsten Tepper Garcia
Date:	02/07/2019
'''

import config.phys_consts as pc
from astropy import units


# Unit definition
MASS = 1. * units.M_sun
LENGTH = 1. * units.kpc
VELOCITY = (1.*units.km/units.s)
TIME = (1.*units.kpc / (units.km / units.s)).decompose().to(units.Gyr)

# Unit name
TIME_UNIT_STR = str(TIME.unit)
MASS_UNIT_STR = str(MASS.unit)
LENGTH_UNIT_STR =  str(LENGTH.unit)
VEL_UNIT_STR = str(VELOCITY.unit).replace(' / ','/')
ANG_MOM_UNIT_STR = str(LENGTH.unit*VELOCITY.unit).replace(' / ','/')
ENERGY_UNIT_STR = str(VELOCITY.unit**2).replace(' / ','/')


# Informative output
# (consider moving this into the units.py module)
print("\nConstants and units:\n")
print("{:>25}{:12.3E} {:15}".format("Gravitational constant:",pc.Grav,"kpc km^2 / Msun s^2"))
print("{:>25}{:12.3f} = {:6E}".format("Mass unit:",MASS,MASS.cgs))
print("{:>25}{:12.3f} = {:6E}".format("Length unit:",LENGTH,LENGTH.cgs))
print("{:>25}{:12.3f} = {:6E}".format("Velocity unit:",VELOCITY,VELOCITY.cgs))
print("{:>25}{:12.3f} = {:6E}".format("Time unit:",TIME,TIME.cgs))
