#import modules
import matplotlib.pyplot as plt
import numpy as np

"""
this program will as the used for the input of initital condidtion, 
and using that condition, plot out the path determined by those initial conditions
"""
# Constants and Initital condititons
g = 9.81  # [m/s^2]
m = 1  # [kg]
l = 0.4  # pendulum lengths in [m]

## Asks the User to pick the initial Condition
print('pick the angles at which you would want to release the double pendulum from')
theta_01 = float(input("Pick a number x between 0-2 such that, theta1 = xpi: ")) 
theta_02 = float(input("Pick a number x between 0-2 such that, theta1 = xpi: ")) 
theta1_0 = np.pi * theta_01 # takes the user inpu to find the initial angle for bob 1
theta2_0 = np.pi * theta_02 # takes the user inpu to find the initial angle for bob 2
omega1_0 = 0 # intital velocity
omega2_0 = 0 # intital velocity

t_0 = 0 # time starts at 0
t_f = 100 # [s]
N = 100000 #stepsize
h = (t_f - t_0) / N # time step

# put all of the initial condition in the same array
r = np.array([theta1_0, omega1_0, theta2_0, omega2_0], float)

# create an array for time
tpoints = np.arange(t_0, t_f, h)

# Empty list to store theta and omega values
theta1_points = []
theta2_points = []
energy_points = []
omega1_points = []
omega2_points = []


# defining a funtion consisting ODEs
def f(r):
    theta1 = r[0]
    omega1 = r[1]
    theta2 = r[2]
    omega2 = r[3]
    
    ##breaking up the equation into smaller sizes (numerator)
    a1 = (omega1 ** 2) * (np.sin(2 * theta1 - 2 * theta2))
    a2 = 2 * (omega2 ** 2) * (np.sin(theta1 - theta2))
    a3 = (g / l) * (np.sin(theta1 - 2 * theta2) + (3 * np.sin(theta1)))

    ##breaking up the equation into smaller sizes (numerator)
    b1 = 4 * (omega1 ** 2) * (np.sin(theta1 - theta2))
    b2 = (omega2 ** 2) * (np.sin(2 * theta1 - 2 * theta2))
    b3 = 2 * (g / l) * (np.sin(2 * theta1 - theta2) - np.sin(theta2))
    
    C = 3 - np.cos(2 * theta1 - 2 * theta2) #denominator
    
    #add it up
    f_omega1 =  -(a1 + a2 + a3)/C
    f_omega2 =  (b1 + b2 + b3)/C
    
    return np.array([omega1, f_omega1, omega2, f_omega2], float)

# Defining the RK4 function
def RK4(r):
    k1 = h * f(r)
    k2 = h * f(r + 0.5 * k1)
    k3 = h * f(r + 0.5 * k2)
    k4 = h * f(r + k3)
    return (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
## Main Loop
for t in tpoints:
    theta1_points.append(r[0])
    theta2_points.append(r[2])
    omega1_points.append(r[1])
    omega2_points.append(r[3])
    r += RK4(r)

# Plots
plt.plot(tpoints, theta1_points, linewidth = 1)
plt.plot(tpoints, theta2_points, linewidth = 1)
plt.title('Motion of a Double Pendulum over Time')
plt.xlabel('Time (s)')
plt.ylabel('$\Theta$')
plt.show()

