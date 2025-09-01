class Fluid:
    def __init__(self, rho=None, k=None, nPow=None, tau0=None, phi0=None, phiInf=None, Tc=None, Ta=None):

        ## Rheology - Modified SMD (Souza Mendes e Dutra (2004)) + Cure(tauY(t))
        # Input Variables
        self.rho = rho  # Density (kg/m³)
        self.k = k  # Consistency Index
        self.nPow = nPow  # Power-law Index
        self.tau0 = tau0 # Dinamic Yield Stress
        self.phi0 = phi0
        self.phiInf = phiInf
        self.Tc = Tc
        self.Ta = Ta

