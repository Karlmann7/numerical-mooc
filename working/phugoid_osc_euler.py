import numpy as np 
import matplotlib.pyplot as plt 

# Final time and time step
T   = 100.0
dt  = 0.02

# initial conditions
z0  = 100.  # altitude
b0  = 10.   # upward velocity resulting from gust
zt  = 100.
g   = 9.81

u   = np.array([z0, b0])

# Number of time step and time grid
N   = int(T/dt) + 1
t   = np.linspace(0.0, T, N)
# t  = np.arange (0.0, T, dt)

# initialize an array to hold the changing elevation values
z   = np.zeros(N)
z[0]= z0  

# time-loop using Euler's method
for n in range(N):
	u    = u + dt*np.array([u[1], g*(1-u[0]/zt)])
	z[n] = u[0]

# Exact solution
z_exact = b0*(zt/g)**0.5*np.sin((g/zt)**.5*t)+\
            (z0-zt)*np.cos((g/zt)**.5*t)+zt

# Plot the solution
plt.figure(figsize=(10,4))   #set plot size
plt.ylim(40,160)             #y-axis plot limits
plt.tick_params(axis='both', labelsize=14) #increase font size for ticks
plt.xlabel('t', fontsize=14) #x label
plt.ylabel('z', fontsize=14) #y label
plt.plot(t, z, 'k-');
plt.plot(t, z_exact)
plt.legend(['Numerical Solution','Analytical Solution']);
plt.show()