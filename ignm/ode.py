import numpy as np

class IgnitionOde:
    def __init__(self,gas,mesh,constants, Tin, Jin, yin):
        
        self.gas = gas
        self.nsp = gas.n_species
        self.MW = gas.molecular_weights
        self.mesh = mesh
        self.cpa = constants[0]
        self.rhoa = constants[1]
        self.da = constants[2]
        self.ya = constants[3]
        self.Tin = Tin
        self.Jin = Jin
        self.yin = yin

    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """

        # State vector is [Y_1, Y_2, ... Y_K, T]
        T = y[self.nsp*self.mesh.n:]
        yj = np.zeros((self.nsp,self.mesh.n))
        wdot = np.zeros((self.nsp,self.mesh.n))
        cp = np.zeros(self.mesh.n)
        cps = np.zeros(self.mesh.n)
        ds = np.zeros(self.mesh.n)
        #rho = np.zeros(self.mesh.n)
        Tprime = np.zeros(self.mesh.n)
        d = np.zeros(self.mesh.n)
        dh = np.zeros((self.nsp,self.mesh.n))
        dyjdt = np.zeros((self.nsp,self.mesh.n))
            
        for i in range(0,self.mesh.n):
            yj[:,i] = y[self.nsp*i:self.nsp*i+self.nsp]
            self.gas.set_unnormalized_mass_fractions(yj[:,i])
            self.gas.TP = T[i], 101325.0
            if T[i]>600:
                dh[:,i] = self.gas.partial_molar_enthalpies
                wdot[:,i] = self.gas.net_production_rates           
            cp[i] = self.gas.cp
            # Diffusion coefficient
            d[i] = self.gas.thermal_conductivity/(self.rhoa*cp[i])
            # Evaluate heat capacity and diffusion coefficient at side boundaries
            cps[i] = (cp[i]+self.cpa)/2
            ds[i] = d[i]*self.da/(self.da+(d[i]-self.da)/2)
            
        Tprime[0] = (-(np.dot(dh[:,0],wdot[:,0])) - self.Jin(t)*cp[0]*(T[0]-self.Tin(t))/self.mesh.dz\
                        + 4*self.rhoa*ds[0]*cps[0]*(300-T[0])/(self.mesh.dx**2))/(self.rhoa*cp[0])
        # Diffusive fluxes at the side boundaries
        jxw = self.rhoa*ds[0]*(yj[:,0]-self.ya)/(self.mesh.dx/2)            
        jxe = self.rhoa*ds[0]*(self.ya-yj[:,0])/(self.mesh.dx/2) 
        
        dyjdt[:,0] = (wdot[:,0]*self.MW - self.Jin(t)*(yj[:,0]-self.yin(t))/self.mesh.dz\
                        -(jxe-jxw)/self.mesh.dx)/(self.rhoa)
        
        for i in range(1,self.mesh.n):
            Tprime[i] = (-(np.dot(dh[:,i],wdot[:,i])) - self.Jin(t)*cp[i]*(T[i]-T[i-1])/self.mesh.dz\
                            + 4*self.rhoa*ds[i]*cps[i]*(300-T[i])/(self.mesh.dx**2))/(self.rhoa*cp[i])
        
            # Diffusive fluxes at the side boundaries
            jxw = self.rhoa*ds[i]*(yj[:,i]-self.ya)/(self.mesh.dx/2)            
            jxe = self.rhoa*ds[i]*(self.ya-yj[:,i])/(self.mesh.dx/2) 
        
            dyjdt[:,i] = (wdot[:,i]*self.MW - self.Jin(t)*(yj[:,i]-self.yin(t))/self.mesh.dz\
                            -(jxe-jxw)/self.mesh.dx)/(self.rhoa)            
        
        print(t)
        return np.hstack((dyjdt.flatten(order='F'), Tprime))