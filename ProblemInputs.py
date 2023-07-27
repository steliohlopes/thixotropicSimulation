class Inputs:
    def __init__(self):
        ##Mesh
        self.meshPath = "PreProcessing/Whistle/"
        self.meshFile = "whistle"

        ## Rheology - Modified SMD (Souza Mendes e Dutra (2004)) + Cure(tauY(t))
        # Input Variables
        self.rho = 1000  # Density (kg/m³)
        self.tau0 = 6.21358  # Dinamic Yield Stress
        self.eta0 = 0.001  # Viscosity Value for Low shear rates
        self.etaInf = 64.1  # Equilibrium Viscosity(Newtonian Plato: Lowgh shear rates)
        self.k = 0.1  # Consistency Index
        self.nPow = 0.572  # Power-law Index
        self.ts = 663  # Caracteristic viscosity buildup time

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
