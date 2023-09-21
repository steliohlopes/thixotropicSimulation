from dolfin import *
import sys
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
        return k*pow(self.gammaDot(u)+DOLFIN_EPS,nPow-1)
    
    def sigmoid(self,field):
        a=500
        H = 1/(1+exp(-a*field))
        return H
    
    def normalized_fluidity(self,phiLocal,phi0,phiInf):
        return (phiLocal-phi0)/(phiInf-phi0)
    
    def phieq(self,k,nPow,phi0,phiInf,u,p,phiLocal,sigmay=0):
        sigmaDev = self.TT(u,p,(1/phiLocal))- ((Identity(len(u)) * tr(self.TT(u,p,(1/phiLocal)))/3))
        # sigmaDev = (nabla_grad(u) + nabla_grad(u).T)/phiLocal
        sigma = pow( inner(sigmaDev,sigmaDev.T)/2 ,0.5) 
        b = pow(abs(sigma-sigmay)/k,1/nPow)/sigma
        dif = phiInf-phi0
        H = self.sigmoid(sigma-sigmay)
        return phi0 + ((dif*b)/(dif+b))*H
    
    def Tc(self):
        tc = 663
        return tc
    
    def Ta(self,dimensionless_phieq):
        ta = 59.2*(pow((1-dimensionless_phieq),1.1)/pow(dimensionless_phieq,0.4))
        return ta
        
    def S(self,dimensionless_phieq):
        s = (8/(exp(dimensionless_phieq/0.09)-1))+1.2
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

        a031=(
                (
                    (self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf)
                     -
            self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf)
                    )
            /self.Tc()
                )
            * 
            self.sigmoid(
                (self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf)-
                self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf)))
            )*self.m


            
        a032=((1-
                self.sigmoid(
                    (self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf)-
                    self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))
                )
               )*
            self.S(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))
                /(self.Ta(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))*self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))
            *
            pow(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf)-
                self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf),
                    (self.S(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))+1)
                    /self.S(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))
                )
            *
            pow(self.normalized_fluidity(self.f,self.fluid.phi0,self.fluid.phiInf),
                    (self.S(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))-1)
                    /self.S(self.normalized_fluidity(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f),self.fluid.phi0,self.fluid.phiInf))
                )
            )*self.m
        
        a03=(a032)*self.mesh.dx(metadata={'quadrature_degree': 3})

        L03 = 0 

        # Complete Weak Form
        F0 = (a01 + a02+a03) - (L01 + L02 + L03)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
