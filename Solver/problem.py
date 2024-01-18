from dolfin import *
import sys
from ufl_legacy import (real,conditional,tanh)

sys.path.append("..")


class Problem:
    def __init__(self, mesh, fluid, boundaries):
        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries

        ##### Functions
        ## Trial and Test function(s)
        self.dw = TrialFunction(self.mesh.functionSpace)
        (self.v, self.q, self.m) = TestFunctions(self.mesh.functionSpace)
        self.w = Function(self.mesh.functionSpace)

        # Split into Velocity and Pressure
        if self.mesh.Dim == 3:
            (self.u, self.p, self.f) = (
                as_vector((self.w[0], self.w[1], self.w[2])),
                self.w[3],
                self.w[4],
            )
        elif self.mesh.Dim == 2:
            (self.u, self.p, self.f) = (
                as_vector((self.w[0], self.w[1])),
                self.w[2],
                self.w[3],
            )

    # Deformation Tensor
    def DD(self, u):
        # Cartesian
        D = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
        return D

    # Deformation Tensor du = 0 ( dx(0) )
    def DD2(self, u):
        #Full Deformation Tensor
        # D = sym(as_tensor([ [        u[0].dx(0)           , (u[1].dx(0) + u[0].dx(1))*0.5 , (u[0].dx(2) + u[2].dx(0))*0.5 ],
        #                     [(u[1].dx(0) + u[0].dx(1))*0.5,         u[1].dx(1)            , (u[2].dx(1) + u[1].dx(2))*0.5 ],
        #                     [(u[0].dx(2) + u[2].dx(0))*0.5, (u[2].dx(1) + u[1].dx(2))*0.5 ,         u[2].dx(2)] ]))


        if self.mesh.Dim == 3:
            D = sym(as_tensor([ [        0        ,     (u[0].dx(1))*0.5         , (u[0].dx(2))*0.5 ],
                                [(u[0].dx(1))*0.5,         u[1].dx(1)            , (u[2].dx(1) + u[1].dx(2))*0.5 ],
                                [(u[0].dx(2))*0.5, (u[2].dx(1) + u[1].dx(2))*0.5 ,         u[2].dx(2)] ]))
            
        elif self.mesh.Dim == 2:
            D = sym(as_matrix([ [        0        ,     (u[0].dx(1))*0.5      ],
                                [(u[0].dx(1))*0.5,         u[1].dx(1)        ],]))

        return D
    
    # Stress Tensor
    def TT(self, u, p, mu):
        # Cartesian
        T = 2 * mu * self.DD(u) - p * Identity(len(u))
        return T
    
    # Stress Tensor DD2
    def TT2(self, u, p, mu):
        # Cartesian
        T = 2 * mu * self.DD2(u) - p * Identity(len(u))
        return T

    def gammaDot(self, u):
        return pow(2 * inner(self.DD(u), self.DD(u)), 0.5)

    def eta(self, k, nPow, u):
        eps = 1e-6
        return k * pow(self.gammaDot(u) + eps, nPow - 1)
    
    def etaSMD(self, k, nPow, u,eta0,etaInf):
        eps = 1e-6
        tauY=DOLFIN_EPS
        gammaDot = self.gammaDot(u)
        return (1-exp(-eta0*gammaDot/tauY))*(tauY/gammaDot + k * pow(gammaDot + eps, nPow - 1) + etaInf)

    def normalized_fluidity(self, phi, phi0, phiInf):
        return (phi - phi0) / (phiInf - phi0)

    def dimensionless_phieq(self, k, nPow, phi0, phiInf, u,localPhi ,sigmay=0):

        V = FunctionSpace(self.mesh.meshObj, self.mesh.Fel)
        # phi  = Function(V)
        # v = TestFunction(V)

        D = self.DD(u)
        gammadot = pow( tr(D * D)/2, 0.5) + DOLFIN_EPS_LARGE

        stress = gammadot/(localPhi*(phiInf-phi0)+phi0)
        max_iter = 100
        erro = 1e-9
        iter = 0
        res = 1

        PHI_v = project(gammadot/(stress),V)
        while iter<=max_iter and res>=erro:
            iter+=1
            g_s=(1/stress)*pow((abs(stress-sigmay)/k),(1/nPow))/((phiInf-phi0)+(1/stress)*pow((abs(stress-sigmay)/k),(1/nPow)))

            c = 500
            H = 0.5*(1+tanh(c*(stress-sigmay)))

            G = g_s*H - (PHI_v-phi0)/(phiInf-phi0) 

            num_der=1/(gammadot)* pow((abs(gammadot/PHI_v-sigmay)/k),(1/nPow)) - pow((k*nPow*PHI_v),(-1)) *pow((abs(gammadot/PHI_v-sigmay)/k),(1/nPow-1))

            den=(phiInf-phi0)+(1/(stress))*pow((abs(stress-sigmay)/k),(1/nPow))

            g_s_der= (phiInf-phi0)*num_der/(den**2)
            H_der=0.5*c*(1-pow((tanh(c*(stress-sigmay))),2))*(-gammadot/PHI_v**2)
            G_der= g_s_der*H+g_s*H_der-1/(phiInf-phi0)
            PHI_N = -(G)/G_der+PHI_v
            PHI_N = project(PHI_N,V)
            stress=gammadot/(PHI_N)
            res = errornorm(PHI_N,PHI_v)
            PHI_v = PHI_N
            begin(f"PhiEq -- r (abs) = {res} (tol = {erro})")

        return (PHI_N-phi0)/(phiInf-phi0)

    def Tc(self):
        tc = 663
        return 1e-6

    def Ta(self, dimensionless_phieq):
        ta = conditional(lt(dimensionless_phieq,1e-7),1e9,59.2 * (
            pow((1 - dimensionless_phieq), 1.1)
            / (pow(dimensionless_phieq, 0.4))
        ))

        return 1e-6

    def S(self, dimensionless_phieq):
        s = conditional(lt(dimensionless_phieq,1e-3),1e9,(8 / (exp(dimensionless_phieq / 0.09) - 1)) + 1.2)
        return s

    
    def GNFEquation(self,model ,wini=None):
        if wini != None:
            self.w = wini

        if model=='newtonian':
            eta = self.fluid.k
        if model=='powerlaw':
            eta = self.eta(self.fluid.k, self.fluid.nPow, self.u)
        if model=='SMD':
            eta = self.etaSMD(self.fluid.k, self.fluid.nPow, self.u,1/self.fluid.phi0,1/self.fluid.phiInf)

        a01 = (
            inner(
                self.TT(self.u, self.p, eta),
                grad(self.v),
            )
        ) * self.mesh.dx()

        outletBCsIndex = tuple(
            self.mesh.subdomains[key]
            for key in self.boundaries.outletBCs
            if key in self.mesh.subdomains
        )
        # Outlet Pressure
        L01 = inner(dot(self.mesh.n , self.TT2(self.u, self.boundaries.Pout , eta)), self.v) * self.mesh.ds(outletBCsIndex)
        
        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(
                self.mesh.subdomains[key]
                for key in self.boundaries.inletBCs
                if key in self.mesh.subdomains
            )
            # Inlet Pressure
            L01 = L01+  inner(dot(self.mesh.n , self.TT2(self.u, self.boundaries.Pin , eta)), self.v) * self.mesh.ds(inletBCsIndex)

        # Mass Conservation(Continuity)
        a02 = (self.q * div(self.u)) * self.mesh.dx()
        L02 = 0

        # Fluidity
        a03 = (
            (eta * (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0) - 1)
            * self.m
            * self.mesh.dx()
        )
        L03 = 0
        # Complete Weak Form
        F0 = (a01 + a02 + a03) - (L01 + L02 + L03)
        # Jacobian Matrix
        J0 = derivative(F0, self.w, self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(
            F0, self.w, self.boundaries.bcs, J0
        )

        return self.problemU0

    def ThixotropicEquation(self, wini=None):
        if wini != None:
            self.w = wini

        a01 = (
            inner(
                self.TT(self.u, self.p, (1 / (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0 ))),
                grad(self.v),
            )
        ) * self.mesh.dx()

        outletBCsIndex = tuple(
            self.mesh.subdomains[key]
            for key in self.boundaries.outletBCs
            if key in self.mesh.subdomains
        )
        # Outlet Pressure
        L01 = inner(dot(self.mesh.n , self.TT2(self.u, self.boundaries.Pout , (1 / (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0 )))), self.v) * self.mesh.ds(outletBCsIndex)
        
        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(
                self.mesh.subdomains[key]
                for key in self.boundaries.inletBCs
                if key in self.mesh.subdomains
            )
            # Inlet Pressure
            L01 = L01+  inner(dot(self.mesh.n , self.TT2(self.u, self.boundaries.Pin , (1 / (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0 )))), self.v) * self.mesh.ds(inletBCsIndex)

        # Mass Conservation(Continuity)
        a02 = (self.q * div(self.u)) * self.mesh.dx()
        L02 = 0

        # Fluidity
        dimensionless_phieq = self.dimensionless_phieq(
                    self.fluid.k,
                    self.fluid.nPow,
                    self.fluid.phi0,
                    self.fluid.phiInf,
                    self.u,
                    self.f
                )
        
        Ta = self.fluid.Ta
        # Ta = self.Ta(dimensionless_phieq)
        Tc = self.fluid.Tc
        # Tc = self.Tc()

        a031 = conditional(
            le(
                self.f,
                dimensionless_phieq,
            ),
            (
                (
                    self.S(
                        dimensionless_phieq
                    )
                    / (
                        Ta
                        * dimensionless_phieq
                    )
                )
                * (
                    pow(
                        (
                                dimensionless_phieq
                                - self.f
                        ),
                        (
                            self.S(
                                dimensionless_phieq
                            )
                            + 1
                        )
                        / self.S(
                            dimensionless_phieq
                        ),
                    )
                )
                * (
                    pow(
                        self.f
                        ,
                        (
                            self.S(
                                dimensionless_phieq
                            )
                            - 1
                        )
                        / self.S(
                            dimensionless_phieq
                        ),
                    )
                )
            ),
            (
                -(
                    self.f
                    - dimensionless_phieq
                )
                / Tc
            ),
        )


        a033 = (
            inner(
                self.u,
                grad(
                    self.f
                ),
            )
        ) 

        a03 = a031 * self.m * self.mesh.dx()
        # L03 = a033 * self.m * self.mesh.dx()

        L03 = dot(self.u*self.m*self.f,self.mesh.n)* self.mesh.ds(outletBCsIndex) - (div(self.u*self.m)*self.f)*self.mesh.dx()

        # Complete Weak Form
        F0 = (a01 + a02 + a03) - (L01 + L02 + L03)

        #SUPG
        #Residual Strong form
        r3 = a031 - a033
        
        fnorm = sqrt(dot(self.f, self.f))
        delta = self.mesh.h/(2.0*fnorm)
        F0 +=delta*r3 *dot(self.u, grad(self.m))*self.mesh.dx()

        # Jacobian Matrix
        J0 = derivative(F0, self.w, self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(
            F0, self.w, self.boundaries.bcs, J0
        )

        return self.problemU0