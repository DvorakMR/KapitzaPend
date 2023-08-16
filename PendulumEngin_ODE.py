import numpy as np
from numpy import sin, cos, pi
from scipy.integrate import odeint
from matplotlib import pyplot as plt
#import matplotlib                  Remove to export to pgf
#matplotlib.use("pgf")
#matplotlib.rcParams.update({
#    "pgf.texsystem": "pdflatex",
#    'font.family': 'serif',
#    'text.usetex': True,
#    'pgf.rcfonts': False,
#})

def find_omega(omegas, y_max):  #Function to find omega for required max y_max
    for omega in omegas:
        phiOpt = odeint(equations, [phi0, x0], time, args=(omega,)) #Diff function to find Optimised phi
        max_value = np.max(phiOpt[:,0])                             #Positive max of phiOpt
        min_value = -np.min(phiOpt[:,0])                            #Max for abs(negative values)
        if abs(max_value) < y_max and min_value < y_max:            #Check if max is below y_max
            return omega
        else:
            print("omega: " , omega, " too small")          #Prints too low omegas for debugging input parameters
#define equations
def equations(y0, t, omega):                    #boundary y0 = [intial angle, angle velocity]
    phi, x = y0
    f = [x, -((R*omega**2 * sin(omega*t)-g)*phi + R*omega**2 * cos(omega*t))/l] #[initial velocity, ode]
    return f

def plot_results(time, phi1):                       #Plotting radian displacement/time
    plt.plot(time, phi1[:,0])
    plt.xlabel("Time (s)")
    plt.ylabel("Angular displacement (rad)")
    #plt.savefig('verysmol.pgf')            #remove to export to pgf
    plt.show()


#parameters
g = 9.81        #m/s^2
rho = 7800      #kg/m^3
L = 1           #Length of pendulum (m)
r1 = 0.04       #radius of sphere (m)
r2 = 0.05       #radius of spherical shell (m)
R = 0.05        #Engine radius (m)
m1 = (4/3)*pi*rho*r1**3             #mass off sphere (kg)
m2 = (4/3)*pi*rho*r2**3 - m1        #mass of spherical shell (kg)    
m = m1+m2                           #total mass of the system
I2 = (2/3)*m2*r2**2                 #Mass moment of inertia of spherical shell
l = (I2 + m*(L+r2)**2)/(m*(L+r2))   #Constant for differential eq.

#inital conditions
initial_angle = 1.0                   #Initial angle speed for pendulum arm (degrees/s)
phi0 = np.radians(initial_angle)      #Degree to radian convertion
initial_angle_speed = 0.0                          #Initial displacement angle for pendulum arm (degrees)
x0 = np.radians(initial_angle_speed)              #Degree to radian convertion
t_max = 100                                 #Max time for plotting on x-axis (s)
time = np.arange(0,t_max, 0.025)            #Step function of time (s)

#find the sol.
omega = 0                                                   #Intialize omega to something 
omegas = np.arange(91,100,1)                           #Array for with stepsize for different omegas
y_max = 0.12                                                   #Maximum allowed displacement for pendulum arm
omega = find_omega(omegas, y_max)                           #Calling find_omega to find required omega for y_max
phi1 = odeint(equations, [phi0, x0], time, args=(omega,))   #Solve differential
max_value = np.max(phi1[:,0])

print(f"Y max below {y_max}: {max_value}, for omega: {omega}") #Print solution parameters
plot_results(time, phi1)    #Function to plot the results
