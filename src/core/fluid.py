class Fluid:
    def __init__(self, config):
        self.config = config.fluid

        # Basic properties
        self.rho = self.config.density
        self.type = self.config.type.value

        # Rheology parameters with sensible fallbacks
        self.k = self.config.k or self.config.viscosity or 1.0
        self.nPow = self.config.n or 1.0
        self.tau0 = getattr(self.config, "tau0", 0.0) or 0.0

        # Structural parameters used by legacy formulations
        self.phi0 = self.config.phi0 or 1e-3
        self.phiInf = self.config.phiInf or max(self.phi0 * 10, self.phi0 + 1e-3)

        # Characteristic times
        self.Ta = self.config.ta or 1.0
        self.Tc = self.config.tc or 1.0

    def __repr__(self) -> str:
        return (
            "Fluid(rho={rho}, k={k}, nPow={n}, tau0={tau0}, "
            "phi0={phi0}, phiInf={phiInf}, Ta={Ta}, Tc={Tc})"
        ).format(
            rho=self.rho,
            k=self.k,
            n=self.nPow,
            tau0=self.tau0,
            phi0=self.phi0,
            phiInf=self.phiInf,
            Ta=self.Ta,
            Tc=self.Tc,
        )