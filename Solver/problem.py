from dolfin import *
import sys
from ufl import tanh
sys.path.append("..")

class Problem:
    def __init__(self, mesh,fluid, boundaries):
        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries

        ##### Functions
        ## Trial and Test function(s)
        self.dw = TrialFunction(self.mesh.functionSpace)
        (self.v, self.q,self.m) = TestFunctions(self.mesh.functionSpace)
        self.w = Function(self.mesh.functionSpace)

        # Split into Velocity and Pressure
        if  self.mesh.Dim == 3:
            (self.u, self.p, self.f) = (as_vector((self.w[0], self.w[1], self.w[2])), self.w[3], self.w[4])
        elif  self.mesh.Dim == 2:
            (self.u, self.p, self.f) = (as_vector((self.w[0], self.w[1])), self.w[2], self.w[3])
    
    # Deformation Tensor
    def DD(self, u):
        # Cartesian
        D = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
        return D


    # Stress Tensor
    def TT(self, u, p, mu):
        # Cartesian
        T = 2 * mu * self.DD(u) - p * Identity(len(u))
        return T


    def gammaDot(self, u):
        return pow(2 * inner(self.DD(u), self.DD(u)), 0.5)


    def eta(self, k,nPow,u):
        eps = 1e-6
        return k*pow(self.gammaDot(u)+eps,nPow-1)
    
    def sigmoid(self,field):
        a=500
        H = 1/(1+exp(-a*field))
        return H
    
    def dsigmoid(self,field):
        a=500
        dH = a*exp(-a*(field))/(pow(exp(-a*field+1),2))
        return dH
    
    def normalized_fluidity(self,phi,phi0,phiInf):
        return (phi-phi0)/(phiInf-phi0)
    
    def func(self,sigma,gammadot,k,nPow,phi0,phiInf,sigmay=0):
        
        b = pow(abs(sigma-sigmay)/k,1/nPow)/(sigma)
        dif = phiInf-phi0
        phieq = b/(dif+b)
        H = self.sigmoid(sigma-sigmay)
        phiV = self.normalized_fluidity(gammadot/sigma,phi0,phiInf)
        return phieq*H-phiV
    
    def dfunc(self,sigma,gammadot,k,nPow,phi0,phiInf,sigmay=0):
        b = pow(abs(sigma-sigmay)/k,1/nPow)/(sigma)
        dif = phiInf-phi0
        phieq = b/(dif+b)
        H = self.sigmoid(sigma-sigmay)

        num = pow(1/k,1/nPow)*(phi0-phiInf)*pow(abs(sigma-sigmay),1/(nPow-2))*(nPow*pow(abs(sigma-sigmay),2)+sigma*(sigmay-sigma))
        den = pow(1/k,1/nPow)*pow(abs(sigma-sigmay),1/(nPow))+sigma*(phiInf-phi0)
        dphieq = num/(nPow*pow(den,2))
        dH = self.dsigmoid(sigma-sigmay)
        dphiV = gammadot/(pow(sigma,2)*(phi0-phiInf))
        return phieq*dH+dphieq*H-dphiV
    
    def dimensionless_phieq(self,k,nPow,phi0,phiInf,u,p,phiLocal,sigmay=0):
        D=self.DD(u)
        gammadot = pow(2*tr(D*D),0.5)+DOLFIN_EPS_LARGE
        sigma = (sigmay + k*pow(gammadot,nPow))+DOLFIN_EPS_LARGE
        # phiV = gammadot/sigma
        # G = self.func(sigma,gammadot,k,nPow,phi0,phiInf,sigmay)
        # dG = self.dfunc(sigma,gammadot,k,nPow,phi0,phiInf,sigmay)
        # phiN = phiV - G/dG

        nIter=0
        maxIter=20
        tol = 1e-9
        res = 2*tol
        V = FunctionSpace(self.mesh.meshObj,self.mesh.Fel)

        while nIter<maxIter and res>tol:
            nIter +=1
            # # phiV = gammadot/sigma
            # phiV = project(gammadot/sigma, V)
            # G = self.func(sigma,gammadot,k,nPow,phi0,phiInf,sigmay)
            # dG = self.dfunc(sigma,gammadot,k,nPow,phi0,phiInf,sigmay)
            # # phiN = phiV - G/dG
            # phiN = project(phiV - G/dG, V)
            # sigma = gammadot/phiN

            # # res = abs(phiN-phiV)
            # res = errornorm(phiN, phiV,'L2')
            # begin(str(res))
            # begin(str(nIter<maxIter and res>tol))

            g_s=(1/sigma)*pow((abs(sigma-sigmay)/k),(1/nPow))/((phiInf-phi0)+(1/sigma)*pow((abs(sigma-sigmay)/k),(1/nPow)))
            c=500           
            H=0.5*(1+tanh(c*(sigma-sigmay)))
            PHI_v=gammadot/(sigma)
            PHI_v = project(PHI_v, V)
            G = g_s*H - (PHI_v-phi0)/(phiInf-phi0)
            num_der=1/(gammadot)*(pow( (abs(gammadot/PHI_v-sigmay)/k),(1/nPow))) - (pow((k*nPow*PHI_v),-1))*(pow((abs(gammadot/PHI_v-sigmay)/k),(1/nPow-1)))
            den=(phiInf-phi0)+(1/(sigma))*pow((abs(sigma-sigmay)/k),(1/nPow))
            g_s_der= (phiInf-phi0)*num_der/(den**2)
            H_der=0.5*c*(1-  pow((tanh(c*(sigma-sigmay))),2) )*(-gammadot/PHI_v**2)
            G_der= g_s_der*H+g_s*H_der-1/(phiInf-phi0)
            PHI_N = -(G)/G_der+PHI_v
            PHI_N = conditional(lt(PHI_N,phi0),phi0,PHI_N)
            PHI_N = conditional(gt(PHI_N,phiInf),phiInf,PHI_N)
            PHI_N = project(PHI_N, V)
            sigma=gammadot/(PHI_N)
            res = errornorm(PHI_N, PHI_v,'L2')

            comm = MPI.comm_world
            if comm.rank == 0:
                begin(f'iteration {nIter}: r (abs) = {res:.2e} (tol = {tol})')
        
        
        return project((PHI_N-phi0)/(phiInf-phi0),V)
    
    def Tc(self):
        tc = 663
        return tc
    
    def Ta(self,dimensionless_phieq):
        ta = 59.2*(pow((1-dimensionless_phieq),1.1)/(pow(dimensionless_phieq+DOLFIN_EPS_LARGE,0.4)))
        return ta
        
    def S(self,dimensionless_phieq):
        s = (8/((exp((dimensionless_phieq+DOLFIN_EPS_LARGE)/0.09)-1)))+1.2
        return s


    def NewtonianEquation(self,wini = None): 
        if wini != None:
            self.w = wini

        a01 = (inner(self.TT(self.u,self.p,self.eta(self.fluid.k,1,self.u)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 

        outletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.outletBCs if key in self.mesh.subdomains)    
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(outletBCsIndex) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.inletBCs if key in self.mesh.subdomains)
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(inletBCsIndex) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        # Fluidity
        a03 = (self.eta(self.fluid.k,1,self.u)*self.f-1)*self.m*self.mesh.dx()
        L03 = 0

        # Complete Weak Form
        F0 = (a01 + a02+a03) - (L01 + L02+L03)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
    
    def PowerLawEquation(self,wini = None):
        if wini != None:
            self.w = wini
        
        a01 = (inner(self.TT(self.u,self.p,self.eta(self.fluid.k,self.fluid.nPow,self.u)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 

        outletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.outletBCs if key in self.mesh.subdomains)                             
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(outletBCsIndex) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.inletBCs if key in self.mesh.subdomains)
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(inletBCsIndex) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0
        
        # Fluidity
        a03 = (self.eta(self.fluid.k,self.fluid.nPow,self.u)*self.f-1)  *self.m*self.mesh.dx()
        L03 = 0
        # Complete Weak Form
        F0 = (a01 + a02+a03) - (L01 + L02+L03)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
    
    def ThixotropicEquation(self,wini = None):
        if wini != None:
            self.w = wini

        a01 = (inner(self.TT(self.u,self.p,(1/self.f)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 

        outletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.outletBCs if key in self.mesh.subdomains)                             
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(outletBCsIndex) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.inletBCs if key in self.mesh.subdomains)
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(inletBCsIndex) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        #Fluidity
        V = FunctionSpace(self.mesh.meshObj,self.mesh.Fel)
        # dimensionless_phieq = self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)
        # normalized_fluidity = project(self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf),V)
        # sigmoid = project(self.sigmoid(normalized_fluidity-dimensionless_phieq),V)
        # S = project(self.S(dimensionless_phieq),V)
        # Ta = project(self.Ta(dimensionless_phieq),V)
        a031=(
                self.sigmoid(self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf)
                             -
                             self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))
                *
                -   (self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf)
                    -
                    self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)
                    ) /self.Tc()
            )*self.m
            
        a032=(
                (1-
                 self.sigmoid(self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf)
                             -
                             self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))
                )
                *
                (
                    self.S(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))
                    /
                    (
                        self.Ta(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))
                        *
                        self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)
                    )
                )
                *
                (
                    pow((abs(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)
                         -
                         self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf))),
                            (self.S(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))+1)
                            /
                           self.S(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))
                        )
                )
                *
                (
                    pow(self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf),
                            (self.S(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))-1)
                            /
                            self.S(self.dimensionless_phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))
                        )
                )
            )*self.m
        
        a033 = (inner(self.u,grad(self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf))))*self.m
        
        a03=(a031+a032+a033)*self.mesh.dx(metadata={'quadrature_degree': 2})


        L03 = 0 

        # Complete Weak Form
        F0 = (a01 + a02+a03) - (L01 + L02 + L03)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
