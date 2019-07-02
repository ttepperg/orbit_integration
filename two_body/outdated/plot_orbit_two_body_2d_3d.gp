#!/opt/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 6    last modified 2019-01-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2018
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')

#----------------------------------------------------------------------------------------

# Author:	Thorsten Tepper Garcia
# Date:	26/06/2019

#----------------------------------------------------------------------------------------
# Invoke example

# Plot 2D data

# $> gnuplot -e 'ndim=2; dataFile="./output/two_body_orbit_int_kepler_2d_leap_phys.dat"; plotRelOrbit = "T"' plot_orbit_two_body.gp

# Plot 3D data in 3D view
# $> gnuplot -e 'ndim=3; dataFile="./output/two_body_orbit_int_kepler_3d_leap_phys.dat"; plotRelOrbit = "T"' plot_orbit_two_body.gp

# Plot 3D data projected onto 2D
# $> gnuplot e 'ndim=3; dataFile="./output/two_body_orbit_int_kepler_3d_leap_phys.dat"; plotRelOrbit = "F"; projPlane = "xz"' plot_orbit_two_body.gp

#----------------------------------------------------------------------------------------
# User settings

# The following are set by the orbit integration program
lengthUnitName = "kpc"
lengthUnit = 1.
velUnitName = "kms^{-1}"
velUnit = 1.
timeUnitName = "Gyr"
timeUnit = 0.978

# EXPERIMENTAL: Analytic solution parameters
plot_ana = 0 # -> 0/1 to turn off/on
semi_latus = 9.1236
ecc = 0.
orbit_angle = 153.9966	# in degree

# defaults

if(!exists('ndim')){
	ndim = 3
}

if(!exists('dataFile')){
	print "\nInput datafile name required (full path).\n"
	quit
}

if(!exists("velVectorScale")){
	velVectorScale = 0.1	# for nice arrows on plot
}
if(!exists("timeStep")){
	timeStep = 0.001
}
if(!exists("timeFreq")){
	timeFreq = 10
}
if(!exists("pauseStep")){
	pauseStep = 0.001
}


# If F, will plot individual orbits instead
if(!exists("plotRelOrbit")){
	plotRelOrbit = "T"
}

#----------------------------------------------------------------------------------------
# Internal settings

# Functions
length(x,y,z) = sqrt(x**2+y**2+z**2)
unit(s) = sprintf(" [%s]",s)
max(a,b) = a>b ? a : b
min(a,b) = a<b ? a : b

set macros

# define indices based on datafile column content
timeCoord = 1
angMom = 2
ePot = 3
eKin = 4

plotColsAngMom = \
	sprintf("($%d):(abs(($%d)/columnhead(%d)-1.))", timeCoord, angMom, angMom)
plotColsEnergy = \
	sprintf("($%d):(abs(($%d+$%d)/(columnhead(%d)+columnhead(%d))-1.))", timeCoord, ePot, eKin, ePot, eKin)


if((ndim == 2)||(exists("projPlane"))){ # ndim = 2 or 3D projection onto 2D

	plotCmdOrbit = "pl"

	timeLabelCoords = "first -0.9*plotRange,0.9*plotRange"

	if((exists("projPlane"))){

		x1Coord = eKin+1
		vx1Coord = x1Coord+1
		y1Coord = x1Coord+2
		vy1Coord = x1Coord+3
		z1Coord = x1Coord+4
		vz1Coord = x1Coord+5
		x2Coord = x1Coord+6
		vx2Coord = x1Coord+7
		y2Coord = x1Coord+8
		vy2Coord = x1Coord+9
		z2Coord = x1Coord+10
		vz2Coord = x1Coord+11


		if(projPlane eq "xy"){ # default setting
			xAxisLabel = "x".unit(lengthUnitName)
			yAxisLabel = "y".unit(lengthUnitName)
			# dummy
			zAxisLabel = ""
		} else {
			if(projPlane eq "xz"){
				aux1Coord = y1Coord
				vaux1Coord = vy1Coord
				y1Coord = z1Coord
				vy1Coord = vz1Coord
				z1Coord = aux1Coord
				vz1Coord = vaux1Coord
				aux2Coord = y2Coord
				vaux2Coord = vy2Coord
				y2Coord = z2Coord
				vy2Coord = vz2Coord
				z2Coord = aux2Coord
				vz2Coord = vaux2Coord
				xAxisLabel = "x".unit(lengthUnitName)
				yAxisLabel = "z".unit(lengthUnitName)
				# dummy
				zAxisLabel = ""
			} else {
				if(projPlane eq "yz"){
					aux1Coord = x1Coord
					vaux1Coord = vx1Coord
					x1Coord = y1Coord
					vx1Coord = vy1Coord
					y1Coord = z1Coord
					vy1Coord = vz1Coord
					z1Coord = aux1Coord
					vz1Coord = vaux1Coord
					aux2Coord = x2Coord
					vaux2Coord = vx2Coord
					x2Coord = y2Coord
					vx2Coord = vy2Coord
					y2Coord = z2Coord
					vy2Coord = vz2Coord
					z2Coord = aux2Coord
					vz2Coord = vaux2Coord
					xAxisLabel = "y".unit(lengthUnitName)
					yAxisLabel = "z".unit(lengthUnitName)
					# dummy
					zAxisLabel = ""
				} else {
					print "unknown projection ", projPlane
					quit
				}
			}
		} # projPlane

	} # if-else 3D projection

	if(ndim == 2){	# there is only one projection in 2D, xy; overwrite any previous settings
		x1Coord = eKin+1
		vx1Coord = x1Coord+1
		y1Coord = x1Coord+2
		vy1Coord = x1Coord+3
		x2Coord = x1Coord+4
		vx2Coord = x1Coord+5
		y2Coord = x1Coord+6
		vy2Coord = x1Coord+7
		# dummy
		z1Coord = 0
		vz1Coord = 0
		z2Coord = 0
		vz2Coord = 0

		xAxisLabel = "x".unit(lengthUnitName)
		yAxisLabel = "y".unit(lengthUnitName)

	}

	statColsX = sprintf("($%d-$%d)", x2Coord, x1Coord)
	statColsY = sprintf("($%d-$%d)", y2Coord, y1Coord)
	# dummy
	statColsZ = sprintf("($%d)", 0)

	plotColsRelOrbit_line = \
		sprintf("($%d-$%d):($0 <= i ? $%d-$%d : 1/0)", x2Coord, x1Coord, y2Coord, y1Coord)
	plotColsRelOrbit_dot = \
		sprintf("($%d-$%d):($0 == i ? $%d-$%d : 1/0)", x2Coord, x1Coord, y2Coord, y1Coord)
	plotColsRelVel = \
		sprintf("($%d-$%d):($0 == i ? $%d-$%d : 1/0):(%f*($%d-$%d)):(%f*($%d-$%d))", \
		x2Coord, x1Coord, y2Coord, y1Coord, \
		velVectorScale, vx2Coord, vx1Coord, velVectorScale, vy2Coord, vy1Coord)
	plotColsRelDist = \
		sprintf("($%d):(length($%d-$%d,$%d-$%d,$%d-$%d))", \
			timeCoord, x2Coord, x1Coord, y2Coord, y1Coord, z2Coord, z1Coord)
	plotColsRelSpeed = \
		sprintf("($%d):(length($%d-$%d,$%d-$%d,$%d-$%d))", \
			timeCoord, vx2Coord, vx1Coord, vy2Coord, vy1Coord, vz2Coord, vz1Coord)

	statColsX1 = sprintf("($%d)", x1Coord)
	statColsY1 = sprintf("($%d)", y1Coord)
	statColsX2 = sprintf("($%d)", x2Coord)
	statColsY2 = sprintf("($%d)", y2Coord)
	# dummy
	statColsZ1 = sprintf("($%d)", 0)
	statColsZ2 = sprintf("($%d)", 0)

	plotColsOrbit1_line = sprintf("($%d):($0 <= i ? $%d : 1/0)", x1Coord, y1Coord)
	plotColsOrbit1_dot = sprintf("($%d):($0 == i ? $%d : 1/0)", x1Coord, y1Coord)
	plotColsVel1 = \
		sprintf("($%d):($0 == i ? $%d : 1/0):(%f*($%d)):(%f*($%d))", \
		x1Coord, y1Coord, velVectorScale, vx1Coord, velVectorScale, vy1Coord)

	plotColsOrbit2_line = sprintf("($%d):($0 <= i ? $%d : 1/0)", x2Coord, y2Coord)
	plotColsOrbit2_dot = sprintf("($%d):($0 == i ? $%d : 1/0)", x2Coord, y2Coord)
	plotColsVel2 = \
		sprintf("($%d):($0 == i ? $%d : 1/0):(%f*($%d)):(%f*($%d))", \
		x2Coord, y2Coord, velVectorScale, vx2Coord, velVectorScale, vy2Coord)


}  # if-else dim = 2 or 3D projection

if((ndim == 3)&&(!exists("projPlane"))){ # full 3D

	plotCmdOrbit = "spl"

	timeLabelCoords = "screen 0.1, 0.9"

	x1Coord = eKin+1
	vx1Coord = x1Coord+1
	y1Coord = x1Coord+2
	vy1Coord = x1Coord+3
	z1Coord = x1Coord+4
	vz1Coord = x1Coord+5
	x2Coord = x1Coord+6
	vx2Coord = x1Coord+7
	y2Coord = x1Coord+8
	vy2Coord = x1Coord+9
	z2Coord = x1Coord+10
	vz2Coord = x1Coord+11

	xAxisLabel = "x".unit(lengthUnitName)
	yAxisLabel = "y".unit(lengthUnitName)
	zAxisLabel = "z".unit(lengthUnitName)

	statColsX = sprintf("($%d-$%d)", x2Coord, x1Coord)
	statColsY = sprintf("($%d-$%d)", y2Coord, y1Coord)
	statColsZ = sprintf("($%d-$%d)", z2Coord, z1Coord)

	plotColsRelOrbit_line = \
		sprintf("($%d-$%d):($0 <= i ? $%d-$%d : 1/0):($%d-$%d)", \
			x2Coord, x1Coord, y2Coord, y1Coord, z2Coord, z1Coord)
	plotColsRelOrbit_dot = \
		sprintf("($%d-$%d):($0 == i ? $%d-$%d : 1/0):($%d-$%d)", \
			x2Coord, x1Coord, y2Coord, y1Coord, z2Coord, z1Coord)
	plotColsRelVel = \
		sprintf("($%d-$%d):($0 == i ? $%d-$%d : 1/0):($%d-$%d):(%f*($%d-$%d)):(%f*($%d-$%d)):(%f*($%d-$%d))", \
			x2Coord, x1Coord, y2Coord, y1Coord, z2Coord, z1Coord, \
			velVectorScale, vx2Coord, vx1Coord, velVectorScale, vy2Coord, vy1Coord, \
			velVectorScale, vz2Coord, vz1Coord)
	plotColsRelDist = \
		sprintf("($%d):(length($%d-$%d,$%d-$%d,$%d-$%d))", \
			timeCoord, x2Coord, x1Coord, y2Coord, y1Coord, z2Coord, z1Coord)
	plotColsRelSpeed = \
		sprintf("($%d):(length($%d-$%d,$%d-$%d,$%d-$%d))", \
			timeCoord, vx2Coord, vx1Coord, vy2Coord, vy1Coord, vz2Coord, vz1Coord)

	statColsX1 = sprintf("($%d)", x1Coord)
	statColsY1 = sprintf("($%d)", y1Coord)
	statColsX2 = sprintf("($%d)", x2Coord)
	statColsY2 = sprintf("($%d)", y2Coord)
	statColsZ1 = sprintf("($%d)", z1Coord)
	statColsZ2 = sprintf("($%d)", z2Coord)

	plotColsOrbit1_line = sprintf("($%d):($0 <= i ? $%d : 1/0):($%d)", x1Coord, y1Coord, z1Coord)
	plotColsOrbit1_dot = sprintf("($%d):($0 == i ? $%d : 1/0):($%d)", x1Coord, y1Coord, z1Coord)
	plotColsVel1 = \
		sprintf("($%d):($0 == i ? $%d : 1/0):($%d):(%f*($%d)):(%f*($%d)):(%f*($%d))", \
		x1Coord, y1Coord, z1Coord, velVectorScale, vx1Coord, velVectorScale, vy1Coord, velVectorScale, vz1Coord)

	plotColsOrbit2_line = sprintf("($%d):($0 <= i ? $%d : 1/0):($%d)", x2Coord, y2Coord, z2Coord)
	plotColsOrbit2_dot = sprintf("($%d):($0 == i ? $%d : 1/0):($%d)", x2Coord, y2Coord, z2Coord)
	plotColsVel2 = \
		sprintf("($%d):($0 == i ? $%d : 1/0):($%d):(%f*($%d)):(%f*($%d)):(%f*($%d))", \
		x2Coord, y2Coord, z2Coord, velVectorScale, vx2Coord, velVectorScale, vy2Coord, velVectorScale, vz2Coord)

} # if-else full 3D


#----------------------------------------------------------------------------------------
# PLOTS

# Energy conservation
set terminal x11 0 persist title "Energy" size 800,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "|{/Symbol D} E(t)/E(0) |"
pl dataFile u @plotColsEnergy w l not


# Angular momentum conservation
set terminal x11 1 persist title "Angular momentum" size 800,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "|{/Symbol D} h(t)/h(0)|"
pl dataFile u @plotColsAngMom w l not

# Orbital history

# Relative distance and relative speed 
set terminal x11 2 persist title "Distance and speed" size 800,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "Distance".unit(lengthUnitName)
set ytics nomirror
set y2label "Speed".unit(velUnitName) rotate by -90
set y2tics nomirror
pl dataFile u @plotColsRelDist w l t "Relative distance", \
'' u @plotColsRelSpeed axes x1y2 w l lw 2 t "Relative speed"
unset y2label
unset y2tics
set ytics mirror


# Orbits

if(plotRelOrbit eq "T"){

	stats dataFile u @statColsX nooutput name 'X'
	stats dataFile u @statColsY nooutput name 'Y'
	Z_min = 0.
	Z_max = 0.
	if((ndim == 3)&&(!exists("projPlane"))){
		stats dataFile u @statColsZ nooutput name 'Z'
	}

	X = max(abs(X_min),abs(X_max))
	Y = max(abs(Y_min),abs(Y_max))
	Z = max(abs(Z_min),abs(Z_max))

	plotRange = 1.1*max(abs(Z),max(abs(X),abs(Y)))

	NumRecords = X_records

} else {

	stats dataFile u @statColsX1 nooutput name 'X1'
	stats dataFile u @statColsY1 nooutput name 'Y1'
	Z1_min = 0.
	Z1_max = 0.
	if((ndim == 3)&&(!exists("projPlane"))){
		stats dataFile u @statColsZ1 nooutput name 'Z1'
	}

	X1 = max(abs(X1_min),abs(X1_max))
	Y1 = max(abs(Y1_min),abs(Y1_max))
	Z1 = max(abs(Z1_min),abs(Z1_max))

	plotRange1 = max(abs(Z1),max(abs(X1),abs(Y1)))

	stats dataFile u @statColsX2 nooutput name 'X2'
	stats dataFile u @statColsY2 nooutput name 'Y2'
	Z2_min = 0.
	Z2_max = 0.
	if((ndim == 3)&&(!exists("projPlane"))){
		stats dataFile u @statColsZ2 nooutput name 'Z2'
	}

	X2 = max(abs(X2_min),abs(X2_max))
	Y2 = max(abs(Y2_min),abs(Y2_max))
	Z2 = max(abs(Z2_min),abs(Z2_max))

	plotRange2 = max(abs(Z2),max(abs(X2),abs(Y2)))

	plotRange = 1.1*max(plotRange1,plotRange2)

	NumRecords = X1_records

} # plot relative orbit?

print "Spatial plot range: ", plotRange

# EXPERIMENTAL: Analytic solution
if(plot_ana){
	set samples 10000
	set terminal x11 3 nopersist
	set polar
	pl semi_latus/(1 + ecc * cos(t - orbit_angle*pi/180.))
	set table "polar.dat"
	repl
	unset table
	unset polar
}

set size square

set xrange [-plotRange:plotRange]
set yrange [-plotRange:plotRange]
if((ndim == 3)&&(!exists("projPlane"))){
	set zrange [-plotRange:plotRange]
}

set xlabel xAxisLabel
set ylabel yAxisLabel
if((ndim == 3)&&(!exists("projPlane"))){
	set zlabel zAxisLabel
}

if(plotRelOrbit eq "T"){

	# Relative orbit
	set terminal x11 3 persist title "Relative Orbit" size 600,600 font "Times-Roman,14" enhanced solid

	if(plot_ana){

		do for [i=NumRecords-2:NumRecords-1]{ \
			unset label
			set label sprintf("T = %5.3f %s", i*timeUnit*timeStep*timeFreq,timeUnitName) \
			at @timeLabelCoords
			@plotCmdOrbit \
			'+' u (0):(0):(0) w p pt 7 ps 2 lw 2 lc rgb "black" t "body 1", \
			dataFile u @plotColsRelOrbit_dot w p pt 7 ps 2 lc rgb "red"  t "body 2", \
			'' u @plotColsRelOrbit_line w l lw 2 t "rel. orbit", \
			'' u @plotColsRelVel w vectors t "rel. vel", \
			'polar.dat' u 1:2 w l lc rgb "black" t "analytic"; \
			pause pauseStep \
		} # do

	} else {

		do for [i=0:NumRecords-1]{ \
			unset label
			set label sprintf("T = %5.3f %s", i*timeUnit*timeStep*timeFreq,timeUnitName) \
			at @timeLabelCoords
			@plotCmdOrbit \
			'+' u (0):(0):(0) w p pt 7 ps 2 lw 2 lc rgb "black" t "body 1", \
			dataFile u @plotColsRelOrbit_dot w p pt 7 ps 2 lc rgb "red"  t "body 2", \
			'' u @plotColsRelOrbit_line w l lw 2 t "rel. orbit", \
			'' u @plotColsRelVel w vectors t "rel. vel"; \
			pause pauseStep \
		} # do

	} # plot analytic solution?

} else {

	# Orbit of each body
	set terminal x11 3 persist title "Individual Orbits" size 600,600 font "Times-Roman,14" enhanced solid

	do for [i=0:NumRecords-1]{ \
		unset label
		set label sprintf("T = %5.3f %s", i*timeUnit*timeStep*timeFreq,timeUnitName) \
		at @timeLabelCoords
		@plotCmdOrbit \
		dataFile u @plotColsOrbit1_dot w p pt 7 ps 2 lc rgb "black" t "body 1", \
		'' u @plotColsOrbit1_line w l t "orbit 1", \
		'' u @plotColsVel1 w vectors t "vel 1", \
		'' u @plotColsOrbit2_dot w p pt 7 ps 2 lc rgb "red" t "body 2", \
		'' u @plotColsOrbit2_line w l t "orbit 2", \
		'' u @plotColsVel2 w vectors t "vel 2"; \
		pause pauseStep \
	} # do

} # plot relative orbit?


# remove temporary file
if(plot_ana){ system("rm -f polar.dat") }

#EOF


