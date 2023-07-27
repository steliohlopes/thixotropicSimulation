class Inputs:
    def __init__(self):
        ##Mesh
        self.meshPath = "PreProcessing/Whistle/"
        self.meshFile = "whistle"

        ## Solver Parameters
        self.nonlinearSolver = "newton"
        self.absTol = 1e-9
        self.relTol = 1e-10
        self.maxIter = 30
        self.linearSolver = "mumps"

        ## No slip Boundaries
        self.noSlipBCs = []
        self.noSlipBCs.append("Wall")

        ##Outlet boundary conditions
        self.Pout = 0
        self.outletBCs = []
        self.outletBCs.append("Outlet")

    def pressureBC(self):
        # Pressure Difference
        self.Pin = 0.1
        self.inletBCs = []
        self.inletBCs.append("Inlet")

    def VelocityBC(self):
        # Constant velocity
        self.Uin = 0.0025
        self.inletBCs = []
        self.inletBCs.append("Inlet")

class Fluid:
    def __init__(self, rho=None, k=None, nPow=None, tau0=None, eta0=None, etaInf=None, ts=None):

        ## Rheology - Modified SMD (Souza Mendes e Dutra (2004)) + Cure(tauY(t))
        # Input Variables
        self.rho = rho  # Density (kg/m³)
        self.k = k  # Consistency Index
        self.nPow = nPow  # Power-law Index
        self.tau0 = tau0 # Dinamic Yield Stress
        self.eta0 = eta0  # Viscosity Value for Low shear rates
        self.etaInf = etaInf  # Equilibrium Viscosity(Newtonian Plato: Lowgh shear rates)
        self.ts = ts  # Caracteristic viscosity buildup time

    @classmethod
    def from_file(cls, file_path):
        #.txt file exemple

        # rho=1000
        # k=0.1
        # nPow=0.572
        # tau0=6.21358
        # eta0=0.001
        # etaInf=64.1
        # ts=663

        with open(file_path, 'r') as file:
            lines = file.readlines()

        attributes = {}
        for line in lines:
            key, value = line.strip().split('=')
            attributes[key.strip()] = float(value.strip())

        return cls(**attributes)
