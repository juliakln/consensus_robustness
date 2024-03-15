# Voter Model
# PySCeS Implementation
# Author: Eric Lofgren (Eric.Lofgren@gmail.com)

Modelname: Voter
Description: PySCes Model Description Language Implementation of SIR model

# Set model to run with numbers of individuals
Species_In_Conc: False	# whether species symbols in rate equations represent concentration (True) or amount (False)
Output_In_Conc: False   # output results of numerical operations in concentrations (True) or amounts (False)

#FIX: Zx Zy

# Differential Equations as Reactions
R1:
	Y > X
	q1*X*Y/(X+Y+20)	# proportion of individuals that come into contact and cause infection
	
R2:
	X > Y
	q2*X*Y/(X+Y+20)

R3:
	Y > X
	noise*Y

R4:
	X > Y
	noise*X

R5:
	Y > X
	q1*Y*Zx/(X+Y+20)

R6:
	X > Y
	q2*X*Zy/(X+Y+20)

# Parameter Values
X = 40
Y = 40
Zx = 10
Zy = 10
q1 = 1
q2 = 1
noise = 0

# Total population size, N
!F N = X+Y+20 	# assignment rule: set value before ODEs are evaluated
