import numpy as np
import cantera as ct
import pandas as pd
import yaml
from scipy.integrate import ode
from numpy import genfromtxt
from scipy.interpolate import interp1d


def solver(input_filename,ts,yin,j0):
    """ Solve the 1-D ignition model """

    class Mesh:     
        def __init__(self,n,z,x):
            
            # Number of grid cells 
            self.n = n   
            self.dz = z/n
            self.dx = x
            # Grid cross-sectional area
            self.a = dx**2 
            self.dv = self.a * self.dz 

    class ignitionOde:
        def __init__(self,gas,mesh):
            
            self.gas = gas
            self.nsp = gas.n_species
            self.MW = gas.molecular_weights
            self.mesh = mesh
            
        def __call__(self, t, y):
            """the ODE function, y' = f(t,y) """

            # State vector is [Y_1, Y_2, ... Y_K, T]
            T = y[self.nsp*self.mesh.n:]
            yj = np.zeros((self.nsp,self.mesh.n))
            wdot = np.zeros((self.nsp,self.mesh.n))
            cp = np.zeros(self.mesh.n)
            cps = np.zeros(self.mesh.n)
            ds = np.zeros(self.mesh.n)
            rho = np.zeros(self.mesh.n)
            Tprime = np.zeros(self.mesh.n)
            d = np.zeros(self.mesh.n)
            dh = np.zeros((self.nsp,self.mesh.n))
            dyjdt = np.zeros((self.nsp,self.mesh.n))
                
            for i in range(0,self.mesh.n):
                yj[:,i] = y[self.nsp*i:self.nsp*i+self.nsp]
                gas.set_unnormalized_mass_fractions(yj[:,i])
                gas.TP = T[i], 101325.0
                if T[i]>600:
                    dh[:,i] = gas.partial_molar_enthalpies
                    wdot[:,i] = gas.net_production_rates           
                cp[i] = gas.cp
                # Diffusion coefficient
                d[i] = gas.thermal_conductivity/(rhoa*cp[i])
                # Evaluate heat capacity and diffusion coefficient at side boundaries
                cps[i] = (cp[i]+cpa)/2
                ds[i] = d[i]*da/(da+(d[i]-da)/2)
                
            Tprime[0] = (-(np.dot(dh[:,0],wdot[:,0])) - Jin(t)*cp[0]*(T[0]-Tin(t))/self.mesh.dz\
                         + 4*rhoa*ds[0]*cps[0]*(300-T[0])/(self.mesh.dx**2))/(rhoa*cp[0])
            # Diffusive fluxes at the side boundaries
            jxw = rhoa*ds[0]*(yj[:,0]-ya)/(self.mesh.dx/2)            
            jxe = rhoa*ds[0]*(ya-yj[:,0])/(self.mesh.dx/2) 
            
            dyjdt[:,0] = (wdot[:,0]*self.MW - Jin(t)*(yj[:,0]-yin(t))/self.mesh.dz\
                          -(jxe-jxw)/self.mesh.dx)/(rhoa)
            
            for i in range(1,self.mesh.n):
                Tprime[i] = (-(np.dot(dh[:,i],wdot[:,i])) - Jin(t)*cp[i]*(T[i]-T[i-1])/self.mesh.dz\
                             + 4*rhoa*ds[i]*cps[i]*(300-T[i])/(self.mesh.dx**2))/(rhoa*cp[i])
            
                # Diffusive fluxes at the side boundaries
                jxw = rhoa*ds[i]*(yj[:,i]-ya)/(self.mesh.dx/2)            
                jxe = rhoa*ds[i]*(ya-yj[:,i])/(self.mesh.dx/2) 
            
                dyjdt[:,i] = (wdot[:,i]*self.MW - Jin(t)*(yj[:,i]-yin(t))/self.mesh.dz\
                              -(jxe-jxw)/self.mesh.dx)/(rhoa)            
            
            print(t)
            return np.hstack((dyjdt.flatten(order='F'), Tprime))
    
    
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

    # Function defining mass fractions at inlet
    def yin(t):
        yinter = np.zeros(nsp)
        for i in range(0,nsp):
            f3 = interp1d(yin1[:,0], yin1[:,i+1], kind='linear')
            yinter[i] = f3(t)
        return yinter
    
    # Initial condition
    gas = ct.Solution(chem)
    gas.TPY = T0, P, mix0
    
    # Store ambient air properties
    ya = gas.Y
    cpa = gas.cp_mass
    da = gas.thermal_conductivity/(gas.density*cpa)
    rhoa = gas.density

    # Set up objects representing the ODE and the solver
    mesh = Mesh(n,z,x)
    yj0 = np.zeros((gas.n_species,mesh.n))
    T1 = np.zeros(mesh.n)+T0
    for i in range(0,mesh.n):
        yj0[:,i] = ya
    y0 = np.hstack((yj0.flatten(order='F'), T1))
    y = pd.DataFrame(columns = [item for item in range(1,y0.size+1)]+["time"])
    y.loc[0] = [y0[item] for item in range(0,y0.size)]+[0]
    ode1 = ignitionOde(gas,mesh)
    solver1 = ode(ode1).set_integrator('dvode', method='bdf', rtol=1e-4, atol=1e-7)
    solver1.set_initial_value(y0, 0)
    dt = dt_max
    count = 1
    
    # Integrate the equations 
    while solver1.successful() and solver1.t < t_sim:
        solver1.integrate(solver1.t + dt)
        y1.loc[count] = [solver1.y[item] for item in range(0,y0.size)]+[solver1.t]
        count += 1
    return y1