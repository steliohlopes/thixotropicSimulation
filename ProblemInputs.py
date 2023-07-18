class Inputs:
    def __init__(self):
        ##Mesh
        self.meshPath = ""
        self.meshFile = ""

        ## Mesh Elements
        # Velocity
        self.velocityElementfamily = "Lagrange"
        self.velocityElementOrder = 2
        # Pressure
        self.pressureElementfamily = "Lagrange"
        self.pressureElementOrder = 1

        ## Rheology - Modified SMD (Souza Mendes e Dutra (2004)) + Cure(tauY(t))
        # Input Variables
        self.rho = 1000  # Density (kg/m³)
        self.tau0 = 6.21358  # Dinamic Yield Stress
        self.eta0 = 0.001  # Viscosity Value for Low shear rates
        self.etaInf = 64.1  # Equilibrium Viscosity(Newtonian Plato: Lowgh shear rates)
        self.k = 0.1  # Consistency Index
        self.nPow = 0.572  # Power-law Index
        self.ts = 663  # Caracteristic viscosity buildup time

        # Solver Parameters
        self.absTol = 1e-9
        self.relTol = 1e-10
        self.maxIter = 30
        self.linearSolver = "mumps"

        ## No slip Boundaries
        self.noSlipBCs = []
        self.noSlipBCs.append("Wall")
