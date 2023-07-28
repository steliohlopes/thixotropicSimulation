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
