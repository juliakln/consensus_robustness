import numpy as np
import gillespy2
import matplotlib.pyplot as plt
import pickle


# Generate x values from 0 to 100
#x_values = list(range(101))
def compute_reaching_time_CI_zealots(z):
    if z==2:
        return 5.95
    elif z == 4:
        return 6.12
    elif z == 6:
        return 6.28
    elif z == 8:
        return 6.46
    elif z == 10:
        return 6.65
    elif z == 12:
        return 6.85
    elif z == 14:
        return 7.06
    elif z == 16:
        return 7.28
    elif z == 18:
        return 7.51
    elif z == 20:
        return 7.75
    elif z == 22:
        return 8.0
    elif z == 24:
        return 8.25
    elif z == 26:
        return 8.51
    elif z == 28: 
        return 8.77
    elif z == 30:
        return 9.02
    elif z == 32:
        return 9.27
    elif z == 34:
        return 9.51
    elif z == 36:
        return 9.74
    elif z == 38:
        return 9.95
    elif z == 40:
        return 10.13
    elif z == 42:
        return 10.28
    elif z == 44:
        return 10.40
    elif z == 50:
        return 10.57
    elif z == 60:
        return 10.48
    elif z == 70:
        return 12.04
    elif z == 80:
        return 27.82
    elif z == 82:
        return 39.94
    elif z == 84:
        return 64.95
    elif z == 86:
        return 128.85
    elif z == 88:
        return 0
    elif z==90:
        return 0
    else:
        return 0

def compute_reaching_time_CI_contrarians(z):
    if z==2:
        return  0 
    elif z == 4:
        return 0
    elif z == 6:
        return 0
    elif z == 8:
        return 0
    elif z == 10:
        return 0
    elif z == 12:
        return 0
    elif z == 14:
        return 8287.98
    elif z == 16:
        return 1004.365
    elif z == 18:
        return 223.17
    elif z == 20:
        return 94.208
    elif z == 22:
        return 52.99
    elif z == 24:
        return 41.49
    elif z == 26:
        return 36.31
    elif z == 28: 
        return 32.02
    elif z == 30:
        return 26.755
    elif z == 32:
        return 25.07
    elif z == 34:
        return 20.95
    elif z == 36:
        return 19.53
    elif z == 38:
        return 16.22
    elif z == 40:
        return 15.488
    elif z == 50:
        return 7.969
    elif z == 60:
        return 3.95
    elif z == 70:
        return 2.333
    elif z == 80:
        return 1.726
    elif z==90:
        return 1.5
    else:
        return 0

def compute_reaching_time_CI_both(z):
    # equall number of both
    if z==2:
        return  0 
    elif z == 4:
        return 0
    elif z == 6:
        return 0
    elif z == 8:
        return 0
    elif z == 10:
        return 0
    elif z == 12:
        return 
    elif z == 14:
        return 0
    elif z == 16:
        return 120509.65
    elif z == 18: # Zx=Zy=5,Cx=Cy=3
        return 0
    elif z == 20:
        return 1612.04
    elif z == 22:
        return 0
    elif z == 24:
        return 157.01
    elif z == 26:
        return 0
    elif z == 28: 
        return 58.51
    elif z == 30:
        return 0
    elif z == 32:
        return 39.798
    elif z == 34:
        return 0
    elif z == 36:
        return 31.9
    elif z == 38:
        return 0
    elif z == 40:
        return 25.96
    elif z == 44:
        return 21.21
    elif z == 50:
        return 0
    elif z == 60:
        return 9.89
    elif z == 64:
        return 7.97
    elif z == 70:
        return 0
    elif z == 80:
        return 4.135
    elif z==90:
        return 0
    else:
        return 0

#x_values = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,50,60,70,80,90]
x_values = list(range(101))
#x_values = list(range(40,101))
#x_values = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Compute y values based on some function or provide your own y values
# For demonstration, let's use a simple function y = x^2
y_values = [compute_reaching_time_CI_zealots(x) for x in x_values]
z_values = [compute_reaching_time_CI_contrarians(x) for x in x_values]
zy_values = [compute_reaching_time_CI_both(x) for x in x_values]

# Plot the values
plt.plot(x_values, y_values, label='CI model with zealots')
plt.plot(x_values, z_values, label='CI model with contrarians')
plt.plot(x_values, zy_values, label='CI model with contrarians and zealots mixed')

# Adding labels and title
plt.xlabel('x values')
plt.ylabel('y values')
plt.title('Plot of expected time to have majority of X with zealots/contrarians varying')
plt.legend()

# Show the plot
plt.savefig('test_reaching.png')
plt.show()

