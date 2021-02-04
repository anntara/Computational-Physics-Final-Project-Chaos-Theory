# import modules
import matplotlib.pyplot as plt
import numpy as np

"""
this program will as the used for the input of initital condidtion, 
and using that condition, plot out the path determined by those initial conditions
"""

# Constants and Initital condititons
G = 1 # gravitational constan, using one to simplify the program
m_A = 150 # [10^30 kg]
m_B = 200 # [10^30 kg]
m_C = 250 # [10^30 kg]

## Asks the User to input the initial Condition for positions
print ('Imagine 3 stars existing on a Eucledian space, input the values in for x and y and assume their distance is in Astronomical Units')
print ('Since 2 stars cannot occupy the same position DO NOT ented the same x and y values for stars A, B and C ')
x_A_0 = float(input("Pick intial position in x direction for A: ")) # initial position on the x axis for  A
y_A_0 = float(input("Pick intial position in y direction for A: ")) # initial position on the y axis for  A
x_B_0 = float(input("Pick intial position in x direction for B: ")) # initial position on the x axis for  B
y_B_0 = float(input("Pick intial position in y direction for B: ")) # initial position on the y axis for  B
x_C_0 = float(input("Pick intial position in x direction for C: ")) # initial position on the x axis for  C
y_C_0 = float(input("Pick intial position in y direction for C: ")) # initial position on the y axis for  C

vx_A_0 = 0 # initial velocity in x direction for A
vy_A_0 = 0 # initial velocity in x direction for A
vx_B_0 = 0 # initial velocity in x direction for A
vy_B_0 = 0 # initial velocity in x direction for A
vx_C_0 = 0 # initial velocity in x direction for A
vy_C_0 = 0 # initial velocity in x direction for A

t_0 = 0 # time starts at 0
t_f = 10  # [Million Years]
N = 400500
h = (t_f - t_0) / N  # initial step size

# put all of the initial condition in the same array
r = np.array([x_A_0, vx_A_0, y_A_0, vy_A_0, 
              x_B_0, vx_B_0, y_B_0, vy_B_0, 
              x_C_0, vx_C_0, y_C_0, vy_C_0], float)

# Empty list to store time and position values
tpoints = []
xpointsA = []
ypointsA = []
xpointsB = []
ypointsB = []
xpointsC = []
ypointsC = []

# defining a funtion consisting 6 first order ODEs
def f(r, t):
    x_A = r[0]
    vx_A = r[1]
    y_A = r[2]
    vy_A = r[3]
    x_B = r[4]
    vx_B = r[5]
    y_B = r[6]
    vy_B = r[7]
    x_C = r[8]
    vx_C = r[9]
    y_C = r[10]
    vy_C = r[11]
    
    # distance between 2 vectors
    rAB = np.sqrt((x_A - x_B) ** 2 + (y_A - y_B) ** 2)
    rBC = np.sqrt((x_C - x_B) ** 2 + (y_C - y_B) ** 2)
    rCA = np.sqrt((x_A - x_C) ** 2 + (y_A - y_C) ** 2)
    
    fx_A = G * m_B * (x_B - x_A) / rAB ** 3 + G * m_C * (x_C - x_A) / rCA ** 3
    fy_A = G * m_B * (y_B - y_A) / rAB ** 3 + G * m_C * (y_C - y_A) / rCA ** 3

    fx_B = G * m_A * (x_A - x_B) / rAB ** 3 + G * m_C * (x_C - x_B) / rBC ** 3
    fy_B = G * m_A * (y_A - y_B) / rAB ** 3 + G * m_C * (y_C - y_B) / rBC ** 3
    
    fx_C = G * m_A * (x_A - x_C) / rCA ** 3 + G * m_B * (x_B - x_C) / rBC ** 3
    fy_C = G * m_A * (y_A - y_C) / rCA ** 3 + G * m_B * (y_B - y_C) / rBC ** 3
    
    return np.array([vx_A, fx_A, vy_A, fy_A,
                     vx_B, fx_B, vy_B, fy_B,
                     vx_C, fx_C, vy_C, fy_C ], float)

# Defining the Adaptive RK4 function
def ATS(r, t, h):
    def RK4(r, t, h):
        k1 = h * f(r, t)
        k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
        k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
        k4 = h * f(r + k3, t + h)
        return (k1 + 2 * k2 + 2 * k3 + k4) / 6 #fouth order Runge-Kutta

    # perform 2 RK steps of step size h
    s_1 = RK4(r, t, h)
    s_2 = RK4(r + s_1, t + h, h)
    r1 = s_1 + s_2

    # perform 1 RK step with step size 2h
    r2 = RK4(r, t, 2 * h)

    # Compute error estimates
    # error for star A
    x1_A = r1[0]
    y1_A = r1[2]
    x2_A = r2[0]
    y2_A = r2[2]
    error_A = np.sqrt((x1_A - x2_A) ** 2 + (y1_A - y2_A) ** 2) / 30 # error in eucledian space

    # error for star B
    x1_B = r1[4]
    y1_B = r1[6]
    x2_B = r2[4]
    y2_B = r2[6]
    error_B = np.sqrt((x1_B - x2_B) ** 2 + (y1_B - y2_B) ** 2) / 30 # error in eucledian space

    # error for star C
    x1_C = r1[8]
    y1_C = r1[10]
    x2_C = r2[8]
    y2_C = r2[10]
    error_C = np.sqrt((x1_C - x2_C) ** 2 + (y1_C - y2_C) ** 2) / 30 # error in eucledian space

    # Use the largest error
    epsilon = max(error_A, error_B, error_C)

    # Finiding Rho using epsilon
    delta = 0.001  # error per unit time
    rho = h * delta / epsilon

    # Calculate factor to multiply h by
    factor = rho ** (1/4)

    # Update h accordingly
    # If target accuracy met, move on to next step
    if  rho >= 1:
        # update t
        t = t + 2 * h

        # Prevent h from getting too big
        if factor > 2:
            h *= 1.5
        else:
            h *= factor

        # Use local extrapolation to better our estimate of the positions
        r1[0]  += (x1_A - x2_A) / 15
        r1[2]  += (y1_A - y2_A) / 15
        r1[4]  += (x1_B - x2_B) / 15
        r1[6]  += (y1_B - y2_B) / 15
        r1[8]  += (x1_C - x2_C) / 15
        r1[10] += (y1_C - y2_C) / 15
        return r1, h, t
    # If target accuracy not met, must redo step with smaller h
    else:
        return ATS(r, t, factor * h)
    
    
## Main Loop
t = t_0
while(t < t_f):
    tpoints.append(t)
    xpointsA.append(r[0])
    ypointsA.append(r[2])
    xpointsB.append(r[4])
    ypointsB.append(r[6])
    xpointsC.append(r[8])
    ypointsC.append(r[10])
    delta_r, h, t = ATS(r, t, h)
    r += delta_r

## Plots
plt.plot(xpointsA, ypointsA, 'g', linewidth = 1, label = 'A')
plt.plot(xpointsB, ypointsB, 'b', linewidth = 1, label = 'B')
plt.plot(xpointsC, ypointsC, 'r', linewidth = 1, label = 'C')
plt.title('Motion of 3 stars over 10 million years')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()

