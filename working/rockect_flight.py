from math import sin, cos, log, ceil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# model parameters:
g   = 9.8      # gravity in m s^{-2}
m_s = 50.      # weight of the rocket shell kg
rho = 1.091    # average air density kg/m^{3}
r   = 0.5      # m^{2}
A   = np.pi*r**2. # maximum cross sectional area of the rocket
v_e = 325.     # exhaust speed m/s^{2}
C_D = 0.15     # drag coefficient --- or D/L if C_L=1
m_p0= 100.     # itial weight of the rocket propellant at t = 0
h0  = 0.
v0  = 0.

def euler_step(u, f, dt, t):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equations.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u, t)


def f(u, t):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    v   = u[1]

    if t < 5.0:
        m_pt = 20.0
        m_p = m_p0 - m_pt*t
    else:
        m_pt = 0.0
        m_p  = 0.0
    
    # print(t,m_p)
    return np.array([  v,
                      -g + 1./(m_s + m_p)*(m_pt*v_e - 0.5*rho*v*np.abs(v)*A*C_D)])



T =  5.0                          # final time
dt = 0.1                           # time increment
N = int(T/dt) + 1                  # number of time-steps
t = np.linspace(0.0, T, N)

# initialize the array containing the solution for each time-step
u = np.empty((N, 2))
u[0] = np.array([h0, v0])# fill 1st element with initial values

# time loop - Euler method
for n in range(1,N):
    u[n] = euler_step(u[n-1], f, dt, t[n])

x = u[:,0]
y = u[:,1]

print(y[-1],t[-1])

# visualization of the path
plt.figure(figsize=(8,6))
plt.grid(True)
plt.xlabel(r'v', fontsize=18)
plt.ylabel(r'h', fontsize=18)
# plt.title('Glider trajectory, flight time = %.2f' % T, fontsize=18)
plt.plot(t,y, 'k-', lw=2);
plt.show()

