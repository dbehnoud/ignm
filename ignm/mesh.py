class Mesh:     
        def __init__(self,n,z,x):
            
            # Number of grid cells 
            self.n = n   
            self.dz = z/n
            self.dx = x
            # Grid cross-sectional area
            self.a = self.dx**2 
            self.dv = self.a * self.dz