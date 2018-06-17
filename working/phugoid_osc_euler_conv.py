import numpy as np 
import matplotlib.pyplot as plt 


def get_error(z, dt):
    """Returns the error relative to analytical solution using L-1 norm.
    
    Parameters
    ----------
    z : array of float
        numerical solution.
    dt : float
        time increment.
        
    Returns
    -------
    err : float
        L_{1} norm of the error with respect to the exact solution.
    """
    N = len(z)
    t = np.linspace(0.0, T, N)
    
    z_exact = b0*(zt/g)**.5*np.sin((g/zt)**.5*t)+\
                (z0-zt)*np.cos((g/zt)**.5*t)+zt
    
    return dt * np.sum(np.abs(z-z_exact))

# Final time
T   = 100.0

# initial conditions
z0  = 100.  # altitude
b0  = 10.   # upward velocity resulting from gust
zt  = 100.
g   = 9.81

# time-increment array
dt_values = np.array([0.1, 0.05, 0.01, 0.005, 0.001])

# array that will contain solution of each grid
z_values = np.empty_like(dt_values, dtype=np.ndarray)

for i, dt in enumerate(dt_values):
    N = int(T/dt)+1    # number of time-steps
    ### discretize the time using numpy.linspace() ###
    t = np.linspace(0.0, T, N)

    # initial conditions
    u = np.array([z0, b0])
    z = np.empty_like(t)
    z[0] = z0
    
    # time loop - Euler method
    for n in range(1,N):
        ### compute next solution using Euler method ###
        u = u + dt*np.array([u[1], g*(1-u[0]/zt)])
        z[n] = u[0]   # store the elevation at time-step n+1
    
    z_values[i] = z.copy()    # store the total elevation calculation grid i


error_values = np.empty_like(dt_values)

for i, dt in enumerate(dt_values):
    ### call the function get_error() ###
    error_values[i] = get_error(z_values[i], dt)

plt.figure(figsize=(10, 6))
plt.tick_params(axis='both', labelsize=14) #increase tick font size
plt.grid(True)                         #turn on grid lines
plt.xlabel('$\Delta t$', fontsize=16)  #x label
plt.ylabel('Error', fontsize=16)       #y label
plt.loglog(dt_values, error_values, 'ko-')  #log-log plot
plt.axis('equal')                      #make axes scale equally;
plt.show()