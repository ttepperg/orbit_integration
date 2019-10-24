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

# Author:			Thorsten Tepper Garcia

#----------------------------------------------------------------------------------------
# Invoke example

# Plot data evolution in 3D view
# $> gnuplot -e 'dataFile="./output/filename.dat"; plotRelOrbit = "T"' plot_orbit_two_body.gp

# Plot data evolution projected onto 2D
# $> gnuplot e 'dataFile="./output/filename.dat"; plotRelOrbit = "T"; projPlane = "xz"' plot_orbit_two_body.gp

# Plot full data (no evolution) projected onto 2D
# $> gnuplot e 'ffw=1; dataFile="./output/filename.dat"; plotRelOrbit = "T"; projPlane = "xz"' plot_orbit_two_body.gp


# The mandatory input parameters and their possible values are:
#
# dataFile - full path to input data file (string)
# delta_t - value of integration time step (float)
# output_freq - value of time step output frequency (int)
#
#
# The optional input parameters and their possible values are:
# 
# plotRelOrbit - T / F (string); T means that the relative orbit is plotted; otherwise the
#				 individual orbit of each body are plotted
# projPlane - xy, xz, yz, or op (string); if set, plot the orbit on an orthogonal Cartesian
# 			  projection, or onto the orbital plane. The latter includes the analytic result.
# ffw - 1 / 0 (int); 0 yields an animation to the orbit's evolution; 1 fast-forwards in time
# 		and plots the full evolution in one step
# pauseStep - (float) pause in seconds between frames when ffw=0
# 
#----------------------------------------------------------------------------------------
# User settings

# The following are set by the orbit integration program
lengthUnitName = "kpc"
lengthUnit = 1.
velUnitName = "kms^{-1}"
velUnit = 1.
timeUnitName = "Gyr"
timeUnit = 0.978


# Sanity checks and defaults

if(!exists('dataFile')){
	print "\nFull path to specify input datafile name, either directly in file or via command-line e.g.:"
	print "gnuplot -e 'dataFile = \"full-path-to-file\"' plot_two_body_orbit.gp\n"
	quit
}
if(!exists("delta_t")){
	print "\nNeed to specify time step, either in file or via command-line e.g.:"
	print "gnuplot -e 'delta_t = 1.e-4' plot_two_body_orbit.gp\n"
	quit
}
if(!exists("output_freq")){
	print "\nNeed to specify time step, either in file or via command-line e.g.:"
	print "gnuplot -e 'output_freq = 10' plot_two_body_orbit.gp\n"
	quit
}
if(!exists("velVectorScale")){
	velVectorScale = 0.1	# for nice arrows on plot
}
if(!exists("pauseStep")){
	pauseStep = 0.001
}
if(!exists("plotRelOrbit")){
	plotRelOrbit = "T"
}
# if 1, fast-forward animation ans plot full orbit evolution in one step
if(!exists("ffw")){
	ffw = 0
}

#----------------------------------------------------------------------------------------
# Internal settings

# Functions
length(x,y,z) = sqrt(x**2+y**2+z**2)
unit(s) = sprintf(" [%s]",s)
max(a,b) = a>b ? a : b
min(a,b) = a<b ? a : b

# Colours
body1Colour = "black"
body2Colour = "red"



set macros

# define indices based on datafile column content
timeCoord = 1
angMomIndex = 2
ePotIndex = 3
eKinIndex = 4

plotColsAngMom = \
	sprintf("($%d):(1.e2*abs(($%d)/columnhead(%d)-1.))", timeCoord, angMomIndex, angMomIndex)
plotColsEnergy = \
	sprintf("($%d):(1.e2*abs(($%d+$%d)/(columnhead(%d)+columnhead(%d))-1.))", timeCoord, ePotIndex, eKinIndex, ePotIndex, eKinIndex)


if(exists("projPlane")){ # orthogonal 3D projection onto 2D

	plotCmdOrbit = "pl"

	timeLabelCoords = "first -0.9*plotRange,0.9*plotRange"

	if((exists("projPlane"))){

		x1Coord = eKinIndex+1
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


		if(projPlane ne "op"){

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
			} # projPlane xy, xz, or xz?
		} # skip settings in case projPlane = "op"

	} # if-else 3D projection

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


} else { # full 3D

	plotCmdOrbit = "spl"

	timeLabelCoords = "screen 0.1, 0.9"

	x1Coord = eKinIndex+1
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


# Settings for system's evolution projected onto orbital plane, i.e. such that the
# relative angular momentum aligns with the positive z-axis.

# if(projPlane eq "op"){

	plotCmdOrbitProj = "pl"

	timeLabelCoordsProj = "first -0.9*plotRangeProj,0.9*plotRangeProj"

	x1CoordProj = eKinIndex+13
	vx1CoordProj = x1CoordProj+1
	y1CoordProj = x1CoordProj+2
	vy1CoordProj = x1CoordProj+3
	x2CoordProj = x1CoordProj+4
	vx2CoordProj = x1CoordProj+5
	y2CoordProj = x1CoordProj+6
	vy2CoordProj = x1CoordProj+7

	xAxisLabelProj = "x_{proj}".unit(lengthUnitName)
	yAxisLabelProj = "y_{proj}".unit(lengthUnitName)

	statColsXProj = sprintf("($%d-$%d)", x2CoordProj, x1CoordProj)
	statColsYProj = sprintf("($%d-$%d)", y2CoordProj, y1CoordProj)

	plotColsRelOrbitProj_line = \
		sprintf("($%d-$%d):($0 <= i ? $%d-$%d : 1/0)", \
			x2CoordProj, x1CoordProj, y2CoordProj, y1CoordProj)
	plotColsRelOrbitProj_dot = \
		sprintf("($%d-$%d):($0 == i ? $%d-$%d : 1/0)", \
			x2CoordProj, x1CoordProj, y2CoordProj, y1CoordProj)
	plotColsRelVelProj = \
		sprintf("($%d-$%d):($0 == i ? $%d-$%d : 1/0):(%f*($%d-$%d)):(%f*($%d-$%d))", \
		x2CoordProj, x1CoordProj, y2CoordProj, y1CoordProj, \
		velVectorScale, vx2CoordProj, vx1CoordProj, velVectorScale, vy2CoordProj, vy1CoordProj)
	plotColsRelDistProj = \
		sprintf("($%d):(length($%d-$%d,$%d-$%d,0.))", \
			timeCoord, x2CoordProj, x1CoordProj, y2CoordProj, y1CoordProj)
	plotColsRelSpeedProj = \
		sprintf("($%d):(length($%d-$%d,$%d-$%d,0.))", \
			timeCoord, vx2CoordProj, vx1CoordProj, vy2CoordProj, vy1CoordProj)

	statColsX1Proj = sprintf("($%d)", x1CoordProj)
	statColsY1Proj = sprintf("($%d)", y1CoordProj)
	statColsX2Proj = sprintf("($%d)", x2CoordProj)
	statColsY2Proj = sprintf("($%d)", y2CoordProj)

	plotColsOrbit1Proj_line = sprintf("($%d):($0 <= i ? $%d : 1/0)", x1CoordProj, y1CoordProj)
	plotColsOrbit1Proj_dot = sprintf("($%d):($0 == i ? $%d : 1/0)", x1CoordProj, y1CoordProj)
	plotColsVel1Proj = \
		sprintf("($%d):($0 == i ? $%d : 1/0):(%f*($%d)):(%f*($%d))", \
		x1CoordProj, y1CoordProj, velVectorScale, vx1CoordProj, velVectorScale, vy1CoordProj)

	plotColsOrbit2Proj_line = sprintf("($%d):($0 <= i ? $%d : 1/0)", x2CoordProj, y2CoordProj)
	plotColsOrbit2Proj_dot = sprintf("($%d):($0 == i ? $%d : 1/0)", x2CoordProj, y2CoordProj)
	plotColsVel2Proj = \
		sprintf("($%d):($0 == i ? $%d : 1/0):(%f*($%d)):(%f*($%d))", \
		x2CoordProj, y2CoordProj, velVectorScale, vx2CoordProj, velVectorScale, vy2CoordProj)

# } # projection plane = orbital plane

#----------------------------------------------------------------------------------------
# PLOTS

#----------------------------------------------------------------------------------------
# Energy conservation
set terminal x11 0 persist title "Energy" size 800,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "|{/Symbol D} E(t)/E(0) | (%)"
set grid
pl dataFile u @plotColsEnergy w l not
unset grid

#----------------------------------------------------------------------------------------
# Angular momentum conservation
set terminal x11 1 persist title "Angular momentum" size 800,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "|{/Symbol D} h(t)/h(0)| (%)"
set grid
pl dataFile u @plotColsAngMom w l not
unset grid

#----------------------------------------------------------------------------------------
# Orbital history (Relative distance and relative speed)
set terminal x11 2 persist title "Distance and speed" size 800,600 font "Times-Roman,14" enhanced solid
set xlabel "Time".unit(timeUnitName)
set ylabel "Distance".unit(lengthUnitName)
set ytics nomirror
set y2label "Speed".unit(velUnitName) rotate by -90
set y2tics nomirror
set grid
pl dataFile u @plotColsRelDist w l t "Relative distance", \
'' u @plotColsRelSpeed axes x1y2 w l lw 2 t "Relative speed"
unset y2label
unset y2tics
unset grid
set ytics mirror

# Uncomment to dump the orbital history to a table
# set table dataFile."_pos_vel_tab"
# pl dataFile u @plotColsRelDist, \
# '' u @plotColsRelSpeed
# unset table


#----------------------------------------------------------------------------------------
# Orbit

#----------------------------------------------------------------------------------------
# 1: Full 3D or 2D orthogonal projection

if(!exists("projPlane")||(projPlane ne "op")){

	set xrange [ * : * ]
	set yrange [ * : * ]

	if(plotRelOrbit eq "T"){

		stats dataFile u @statColsX nooutput name 'X'
		stats dataFile u @statColsY nooutput name 'Y'
		Z_min = 0.
		Z_max = 0.
		if(!exists("projPlane")){
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
		if(!exists("projPlane")){
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
		if(!exists("projPlane")){
			stats dataFile u @statColsZ2 nooutput name 'Z2'
		}

		X2 = max(abs(X2_min),abs(X2_max))
		Y2 = max(abs(Y2_min),abs(Y2_max))
		Z2 = max(abs(Z2_min),abs(Z2_max))

		plotRange2 = max(abs(Z2),max(abs(X2),abs(Y2)))

		plotRange = 1.1*max(plotRange1,plotRange2)

		NumRecords = X1_records

	} # plot relative orbit?

	print "Spatial plot range (intrinsic orbit): ", plotRange


	set size square

	set xrange [-plotRange:plotRange]
	set yrange [-plotRange:plotRange]
	if(!exists("projPlane")){
		set zrange [-plotRange:plotRange]
	}

	set xlabel xAxisLabel
	set ylabel yAxisLabel
	if(!exists("projPlane")){
		set zlabel zAxisLabel
	}

	if(plotRelOrbit eq "T"){

		# Relative orbit (full 3D or 2D orthogonal projection)
		set terminal x11 3 persist title "Relative Orbit" size 600,600 font "Times-Roman,14" enhanced solid

		if(ffw){

			do for [i=NumRecords-2:NumRecords-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoords
				@plotCmdOrbit \
				'+' u (0):(0):(0) w p pt 7 ps 2 lw 2 lc rgb body1Colour t "body 1", \
				dataFile u @plotColsRelOrbit_dot w p pt 7 ps 2 lc rgb "red"  t "body 2", \
				'' u @plotColsRelOrbit_line w l lw 2 t "rel. orbit", \
				'' u @plotColsRelVel w vectors t "rel. vel"; \
				pause pauseStep \
			} # do

		} else {

			do for [i=0:NumRecords-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoords
				@plotCmdOrbit \
				'+' u (0):(0):(0) w p pt 7 ps 2 lw 2 lc rgb body1Colour t "body 1", \
				dataFile u @plotColsRelOrbit_dot w p pt 7 ps 2 lc rgb "red"  t "body 2", \
				'' u @plotColsRelOrbit_line w l lw 2 t "rel. orbit", \
				'' u @plotColsRelVel w vectors t "rel. vel"; \
				pause pauseStep \
			} # do

		} # fast-forward?

	} else {

		# Orbit of each body (full 3D or 2D orthogonal projection)
		set terminal x11 3 persist title "Individual Orbits" size 600,600 font "Times-Roman,14" enhanced solid

		if(ffw){

			do for [i=NumRecords-2:NumRecords-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoords
				@plotCmdOrbit \
				dataFile u @plotColsOrbit1_dot w p pt 7 ps 2 lc rgb body1Colour t "body 1", \
				'' u @plotColsOrbit1_line w l t "orbit 1", \
				'' u @plotColsVel1 w vectors t "vel 1", \
				'' u @plotColsOrbit2_dot w p pt 7 ps 2 lc rgb body2Colour t "body 2", \
				'' u @plotColsOrbit2_line w l t "orbit 2", \
				'' u @plotColsVel2 w vectors t "vel 2"; \
				pause pauseStep \
			} # do

		} else {

			do for [i=0:NumRecords-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoords
				@plotCmdOrbit \
				dataFile u @plotColsOrbit1_dot w p pt 7 ps 2 lc rgb body1Colour t "body 1", \
				'' u @plotColsOrbit1_line w l t "orbit 1", \
				'' u @plotColsVel1 w vectors t "vel 1", \
				'' u @plotColsOrbit2_dot w p pt 7 ps 2 lc rgb body2Colour t "body 2", \
				'' u @plotColsOrbit2_line w l t "orbit 2", \
				'' u @plotColsVel2 w vectors t "vel 2"; \
				pause pauseStep \
			} # do

		} # fast-forward?

	} # plot relative orbit (full 3D or 2D orthogonal projection)?


} else {
#----------------------------------------------------------------------------------------
# 2: System's evolution projected onto orbital plane

	set xrange [ * : * ]
	set yrange [ * : * ]

	if(plotRelOrbit eq "T"){

		stats dataFile u @statColsXProj nooutput name 'XProj'
		stats dataFile u @statColsYProj nooutput name 'YProj'

		XProj = max(abs(XProj_min),abs(XProj_max))
		YProj = max(abs(YProj_min),abs(YProj_max))

		plotRangeProj = 1.1*max(abs(XProj),abs(YProj))

		NumRecordsProj = XProj_records

	} else {

		stats dataFile u @statColsX1Proj nooutput name 'X1Proj'
		stats dataFile u @statColsY1Proj nooutput name 'Y1Proj'

		X1Proj = max(abs(X1Proj_min),abs(X1Proj_max))
		Y1Proj = max(abs(Y1Proj_min),abs(Y1Proj_max))

		plotRange1Proj = max(abs(X1Proj),abs(Y1Proj))

		stats dataFile u @statColsX2Proj nooutput name 'X2Proj'
		stats dataFile u @statColsY2Proj nooutput name 'Y2Proj'

		X2Proj = max(abs(X2Proj_min),abs(X2Proj_max))
		Y2Proj = max(abs(Y2Proj_min),abs(Y2Proj_max))

		plotRange2Proj = max(abs(X2Proj),abs(Y2Proj))

		plotRangeProj = 1.1*max(plotRange1Proj,plotRange2Proj)

		NumRecordsProj = X1Proj_records

	} # plot relative orbit?

	print "Spatial plot range (orbital plane): ", plotRangeProj


	set size square

	set xrange [-plotRangeProj:plotRangeProj]
	set yrange [-plotRangeProj:plotRangeProj]

	set xlabel xAxisLabelProj
	set ylabel yAxisLabelProj

	if(plotRelOrbit eq "T"){

		# Relative orbit (projected onto orbital plane)
		set terminal x11 4 persist title "Relative orbit (orbital plane)" size 600,600 font "Times-Roman,14" enhanced solid

		if(!(ffw)){

			do for [i=0:NumRecordsProj-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoordsProj
				@plotCmdOrbitProj \
				'+' u (0):(0) w p pt 7 ps 2 lw 2 lc rgb body1Colour t "body 1", \
				dataFile u @plotColsRelOrbitProj_dot w p pt 7 ps 2 lc rgb "red"  t "body 2", \
				'' u @plotColsRelOrbitProj_line w l lw 2 t "rel. orbit", \
				'' u @plotColsRelVelProj w vectors t "rel. vel", \
				'' u 25:26 w l lc rgb "black" t "analytic"; \
				pause pauseStep \
			} # do

		}
		if(ffw){

			do for [i=NumRecordsProj-2:NumRecordsProj-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoordsProj
				@plotCmdOrbitProj \
				'+' u (0):(0) w p pt 7 ps 2 lw 2 lc rgb body1Colour t "body 1", \
				dataFile u @plotColsRelOrbitProj_dot w p pt 7 ps 2 lc rgb "red"  t "body 2", \
				'' u @plotColsRelOrbitProj_line w l lw 2 t "rel. orbit", \
				'' u @plotColsRelVelProj w vectors t "rel. vel", \
				'' u 25:26 w l lc rgb "black" t "analytic" \
			} # do

		}

	} else {

		# Orbit of each body (projected onto orbital plane)
		set terminal x11 4 persist title "Individual Orbits (orbital plane)" size 600,600 font "Times-Roman,14" enhanced solid

		if(ffw){

			do for [i=NumRecordsProj-2:NumRecordsProj-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoordsProj
				@plotCmdOrbitProj \
				dataFile u @plotColsOrbit1Proj_dot w p pt 7 ps 2 lc rgb body1Colour t "body 1", \
				'' u @plotColsOrbit1Proj_line w l t "orbit 1", \
				'' u @plotColsVel1Proj w vectors t "vel 1", \
				'' u @plotColsOrbit2Proj_dot w p pt 7 ps 2 lc rgb body2Colour t "body 2", \
				'' u @plotColsOrbit2Proj_line w l t "orbit 2", \
				'' u @plotColsVel2Proj w vectors t "vel 2" \
			} # do

		} else {
	

			do for [i=0:NumRecordsProj-1]{ \
				unset label
				set label sprintf("T = %7.5f %s", i*timeUnit*delta_t*output_freq,timeUnitName) \
				at @timeLabelCoordsProj
				@plotCmdOrbitProj \
				dataFile u @plotColsOrbit1Proj_dot w p pt 7 ps 2 lc rgb body1Colour t "body 1", \
				'' u @plotColsOrbit1Proj_line w l t "orbit 1", \
				'' u @plotColsVel1Proj w vectors t "vel 1", \
				'' u @plotColsOrbit2Proj_dot w p pt 7 ps 2 lc rgb body2Colour t "body 2", \
				'' u @plotColsOrbit2Proj_line w l t "orbit 2", \
				'' u @plotColsVel2Proj w vectors t "vel 2"; \
				pause pauseStep \
			} # do

		} # fast-forward?
	

	} # plot relative orbit?

} # full 3D, orthogonal projection, or projection onto orbital plane?

#EOF


