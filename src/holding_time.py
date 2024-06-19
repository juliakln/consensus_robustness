import numpy as np
import gillespy2
import matplotlib.pyplot as plt
import pickle


# Generate x values from 0 to 100
#x_values = list(range(101))
def compute_holding_time_CI_zealots(z):
    if z==2:
        return     
    elif z == 4:
        return 81853.46
    elif z == 6:
        return 50839.65
    elif z == 8:
        return 35842.81
    elif z == 10:
        return 26867
    elif z == 12:
        return 20862.20
    elif z == 14:
        return 16546.94
    elif z == 16:
        return 13214.05
    elif z == 18:
        return 9797.61
    elif z == 20:
        return 4242.5
    elif z == 22:
        return 740.66
    elif z == 24:
        return 150
    elif z == 26:
        return 43.11
    elif z == 28: 
        return 14.08
    elif z == 30:
        return 5.58
    elif z == 32:
        return 2.76
    elif z == 34:
        return 1.65
    elif z == 36:
        return 1.15
    elif z == 38:
        return 0.886
    elif z == 40:
        return 0.737
    elif z == 50:
        return 0.49222
    elif z == 60:
        return 0.43
    elif z == 70:
        return 0.3473
    elif z == 80:
        return 0.265
    elif z==90:
        return 0.218
    else:
        return 0
    
def compute_holding_time_hack_CI_zealots(z):
    if z==2:
        return     
    elif z == 4:
        return# //81853.46
    elif z == 6:
        return #//50839.65
    elif z == 8:
        return# //35842.81
    elif z == 10:
        return #//26867
    elif z == 12:
        return #//0.582
    elif z == 14:
        return #//16546.94
    elif z == 16:
        return #//13214.05
    elif z == 18:
        return# //9797.61
    elif z == 20:
        return 326.22 # //4242.5 TOOK REALLY LONG TIME TO COMPUTE
    elif z == 22:
        return 32.063 #//740.66
    elif z == 24:
        return 6.36
    elif z == 26:
        return 1.644; #//43.11
    elif z == 28: 
        return 0.518#//14.08
    elif z == 30:
        return 0.35; #//5.58
    elif z == 32:
        return 0.266;#//2.76
    elif z == 34:
        return 0.229;#//1.65
    elif z == 36:
        return 0.222; #//1.15
    elif z == 38:
        return 0.221;#//0.886
    elif z == 40:
        return 0.216; #//0.737
    elif z == 50:
        return #//0.49222
    elif z == 60:
        return #//0.43
    elif z == 70:
        return #//0.3473
    elif z == 80:
        return #//0.265
    elif z==90:
        return #//0.218
    else:
        return 0

def compute_holding_time_hack_CI_contrarians(z):
    # when starting with 55x and 45y
    if z==2:
        return  0 
    elif z == 4:
        return 0
    elif z == 6:
        return 0
    elif z == 8:
        return 0 #0.122?!
    elif z == 10:
        return 240.06
    elif z == 12:
        return 14.855
    elif z == 14:
        return 0.6
    elif z == 16:
        return 0
    elif z == 18:
        return 0
    elif z == 20:
        return 0.142
    elif z == 22:
        return 0.141
    elif z == 24:
        return 0.13
    elif z == 26:
        return 0.13
    elif z == 28: 
        return 0.129
    elif z == 30:
        return 0.128
    elif z == 32:
        return 0.1296
    elif z == 34:
        return 0.129
    elif z == 36:
        return 0.131
    elif z == 38:
        return 0.130
    elif z == 40:
        return 0.128
    elif z == 50:
        return 0.131
    elif z == 60:
        return 0.121
    elif z == 70:
        return 0.121
    elif z == 80:
        return 0.123
    elif z==90:
        return 0.125
    else:
        return 0

def compute_holding_time_hack_CI_both(z):
    # when starting with 55x and 45y
    if z==2:
        return  0 
    elif z == 4:
        return 0
    elif z == 6:
        return 0
    elif z == 8:
        return 0 
    elif z == 12:
        return 2466.24
    elif z == 16:
        return 4.035
    elif z == 20:
        return 0.176
    elif z == 24:
        return 0.22
    elif z == 28: 
        return 0.16
    elif z == 32:
        return 0.14
    elif z == 36:
        return 0
    elif z == 40:
        return 0
    elif z == 48:
        return 0
    elif z == 60:
        return 0
    elif z == 80:
        return 0
    elif z==92:
        return 0
    else:
        return 0

#x_values = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,50,60,70,80,90]
x_values = list(range(101))
#x_values = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Compute y values based on some function or provide your own y values
# For demonstration, let's use a simple function y = x^2
yr_values = [compute_holding_time_CI_zealots(x) for x in x_values]
y_values = [compute_holding_time_hack_CI_zealots(x) for x in x_values]
z_values = [compute_holding_time_hack_CI_contrarians(x) for x in x_values]
yz_values = [compute_holding_time_hack_CI_both(x) for x in x_values]


# Plot the values
#plt.plot(x_values, yr_values, label='CI model with zealots real time')
plt.plot(x_values, y_values, label='CI model with zealots')
plt.plot(x_values, z_values, label='CI model with contrarians')
plt.plot(x_values, yz_values, label='CI model with both')

# Adding labels and title
plt.xlabel('x values')
plt.ylabel('y values')
plt.title('Plot of expected holding time with majority X with zealots/contrarians varying')
plt.legend()

# Show the plot
plt.savefig('test_holding.png')
plt.show()
