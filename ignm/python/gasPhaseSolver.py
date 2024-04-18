import numpy as np
import cantera as ct
import pandas as pd
import yaml
from scipy.integrate import ode
from numpy import genfromtxt
from scipy.interpolate import interp1d
from ignm.mesh import Mesh
from ignm.ode import IgnitionOde

# Function defining mass fractions at inlet
def yin(t, nsp, yin1):
    yinter = np.zeros(nsp)
    for i in range(0,nsp):
        f3 = interp1d(yin1[:,0], yin1[:,i+1], kind='linear')
        yinter[i] = f3(t)
    return yinter

def solve(input_filename,ts,yin,j0):
    """ Solve the 1-D gas-phase combustion """    
    # Input data
    with open(input_filename, 'r') as f:
            inputs = yaml.safe_load(f)
    n = inputs['grid_size']
    z = inputs['domain_size_z']
    x = inputs['domain_size_x']
    chem = inputs['mechanism']
    T0 = inputs['initial_temperature']
    P = inputs['pressure']
    mix0 = inputs['mixture']
    t_sim = inputs['simulation_time']
    dt_max = inputs['max_time_step']
    ts1 = genfromtxt(ts, delimiter=',')
    yin1 = genfromtxt(yin, delimiter=',')
    j01 = genfromtxt(j0, delimiter=',')

    # Temperature and mass flux at inlet
    Tin = interp1d(ts1[:,0], ts1[:,1], kind='linear')
    Jin = interp1d(j01[:,0], j01[:,1], kind='linear')
    
    
    # Initial condition
    gas = ct.Solution(chem)
    gas.TPY = T0, P, mix0
    
    # Store ambient air properties
    ya = gas.Y
    cpa = gas.cp_mass
    da = gas.thermal_conductivity/(gas.density*cpa)
    rhoa = gas.density
    constants = (cpa,rhoa,da,ya)

    # Set up objects representing the ODE and the solver
    mesh = Mesh(n,z,x)
    yj0 = np.zeros((gas.n_species,mesh.n))
    T1 = np.zeros(mesh.n)+T0
    for i in range(0,mesh.n):
        yj0[:,i] = ya
    y0 = np.hstack((yj0.flatten(order='F'), T1))
    y = pd.DataFrame(columns = [item for item in range(1,y0.size+1)]+["time"])
    y.loc[0] = [y0[item] for item in range(0,y0.size)]+[0]
    ode1 = IgnitionOde(gas,mesh,constants, Tin, Jin, yin1)
    solver1 = ode(ode1).set_integrator('dvode', method='bdf', rtol=1e-4, atol=1e-7)
    solver1.set_initial_value(y0, 0)
    dt = dt_max
    count = 1
    
    # Integrate the equations 
    while solver1.successful() and solver1.t < t_sim:
        solver1.integrate(solver1.t + dt)
        y.loc[count] = [solver1.y[item] for item in range(0,y0.size)]+[solver1.t]
        count += 1
    return y