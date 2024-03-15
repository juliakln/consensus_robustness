# Cross-Inhibition Model
# PySCeS Implementation
# Author: Eric Lofgren (Eric.Lofgren@gmail.com)

Modelname: CrossInhibition
Description: PySCes Model Description Language Implementation of SIR model

# Set model to run with numbers of individuals
Species_In_Conc: False	# whether species symbols in rate equations represent concentration (True) or amount (False)
Output_In_Conc: False   # output results of numerical operations in concentrations (True) or amounts (False)

#FIX: Zx Zy

# Differential Equations as Reactions
# cross-inhibition
R1:
	Y > U
	q1*X*Y/(X+Y+U+20)	# proportion of individuals that come into contact and cause infection
	
R2:
	X > U
	q2*X*Y/(X+Y+U+20)

# noise type 2
R3:
	Y > X
	noise*Y

R4:
	X > Y
	noise*X

# zealots
R5:
	Y > U
	q1*Y*Zx/(X+Y+U+20)

R6:
	X > U
	q2*X*Zy/(X+Y+U+20)

R7:
	U > X
	q1*U*Zx/(X+Y+U+20)

R8:
	U > Y
	q2*U*Zy/(X+Y+U+20)

# recruitment
R9:
	U > X
	q1*U*X/(X+Y+U+20)

R10:
	U > Y
	q2*U*Y/(X+Y+U+20)

# noise type 1

R11:
	Y > U
	Y*noise
R12:
	U > Y
	U*noise
R13:
	X > U
	X*noise
R14:
	U > X
	U*noise


# Parameter Values
X = 40
Y = 40
U = 0
Zx = 10
Zy = 10
q1 = 1
q2 = 1
noise = 0

# Total population size, N
!F N = X+Y+U+20 	# assignment rule: set value before ODEs are evaluated
