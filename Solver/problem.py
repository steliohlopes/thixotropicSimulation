from dolfin import *
import sys
from ufl_legacy import (real,conditional,tanh)

sys.path.append("..")


class Problem:
    def __init__(self, mesh, fluid, boundaries,L,U):
        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries
        self.L = L
        self.U = U

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
        D = 0.5*self.U * (nabla_grad(u) + nabla_grad(u).T)
        return D

    def DD_direction(self, u,direction):
        #Full Deformation Tensor
        # D = sym(as_tensor([ [        u[0].dx(0)           , (u[1].dx(0) + u[0].dx(1))*0.5 , (u[0].dx(2) + u[2].dx(0))*0.5 ],
        #                     [(u[1].dx(0) + u[0].dx(1))*0.5,         u[1].dx(1)            , (u[2].dx(1) + u[1].dx(2))*0.5 ],
        #                     [(u[0].dx(2) + u[2].dx(0))*0.5, (u[2].dx(1) + u[1].dx(2))*0.5 ,         u[2].dx(2)] ]))
        if self.mesh.Dim == 3:
            if direction==0:
                D = sym(as_tensor([ [        0        ,     (u[0].dx(1))*self.U*0.5         , (u[0].dx(2))*self.U*0.5 ],
                                [(u[0].dx(1))*self.U*0.5,         u[1].dx(1)*self.U            , (u[2].dx(1) + u[1].dx(2))*self.U*0.5 ],
                                [(u[0].dx(2))*self.U*0.5, (u[2].dx(1) + u[1].dx(2))*self.U*0.5 ,         u[2].dx(2)*self.U] ]))
                
            elif direction==1:
                D = sym(as_tensor([ [        u[0].dx(0)*self.U          , (u[1].dx(0))*self.U*0.5 , (u[0].dx(2) + u[2].dx(0))*self.U*0.5 ],
                            [(u[1].dx(0))*self.U*0.5,         0            , ( u[1].dx(2))*self.U*0.5 ],
                            [(u[0].dx(2) + u[2].dx(0))*self.U*0.5, ( u[1].dx(2))*self.U*0.5 ,         u[2].dx(2)*self.U] ]))
                
            elif direction==2:
                D = sym(as_tensor([ [        u[0].dx(0)*self.U           , (u[1].dx(0) + u[0].dx(1))*self.U*0.5 , ( u[2].dx(0))*self.U*0.5 ],
                                [(u[1].dx(0) + u[0].dx(1))*self.U*0.5,         u[1].dx(1)*self.U            , (u[2].dx(1))*self.U*0.5 ],
                                [(      u[2].dx(0))*self.U*0.5,            (u[2].dx(1))*self.U*0.5          ,        0] ]))
                
        elif self.mesh.Dim == 2:
            if direction==0:
                D = sym(as_matrix([ [        0        ,     (u[0].dx(1))*self.U*0.5      ],
                                [(u[0].dx(1))*self.U*0.5,         u[1].dx(1)*self.U        ],]))
                
            elif direction==1:
                D = sym(as_matrix([ [u[0].dx(0)*self.U      , (u[1].dx(0))*self.U*0.5 ],
                                [(u[1].dx(0))*self.U*0.5,         0        ] ]))
                
        return D
    
    # Stress Tensor
    def TT(self, u, p, phi):
        pinf = self.boundaries.Pin
        p0 = self.boundaries.Pout
        dp=pinf-p0
        
        phiInf = self.fluid.phiInf
        phi0 = self.fluid.phi0
        dphi = phiInf - phi0
        # Cartesian
        T = 2/dp * 1/(phi*dphi + phi0) * self.DD(u) - ((p*dp + p0)/dp) * Identity(len(u))
        return T
    
    def TT_direction(self, u, p, phi,direction):
        pinf = self.boundaries.Pin
        p0 = self.boundaries.Pout
        dp=pinf-p0
        
        phiInf = self.fluid.phiInf
        phi0 = self.fluid.phi0
        dphi = phiInf - phi0
        # Cartesian
        T = 2/dp * 1/(phi*dphi + phi0) * self.DD_direction(u,direction) - ((p*dp + p0)/dp) * Identity(len(u))
        
        return T

    #Dimensional GammaDot
    def gammaDot(self, u):
        DD = 0.5 * (nabla_grad(u*self.U) + nabla_grad(u*self.U).T)
        return pow(2 * inner(DD, DD), 0.5)
    
    def etaSMD(self, k, nPow, u,eta0,etaInf):
        eps = 1e-6
        tauY=DOLFIN_EPS
        gammaDot = self.gammaDot(u)
        return (1-exp(-eta0*gammaDot/tauY))*(tauY/gammaDot + k * pow(gammaDot + eps, nPow - 1) + etaInf)

    def normalized_fluidity(self, phi, phi0, phiInf):
        return (phi - phi0) / (phiInf - phi0)

    def dimensionless_phieq(self, k, nPow, phi0, phiInf, u,localPhi ,sigmay=0):
        V = FunctionSpace(self.mesh.meshObj, self.mesh.Fel)
        phi  = Function(V)
        v = TestFunction(V)
        # u = u*self.U
        
        # if self.mesh.Dim == 3:
        #     TraceD2 =  pow(u[0].dx(0),2) + pow(u[1].dx(1),2) + pow( u[2].dx(2),2) + 2*(pow((u[1].dx(0) + u[0].dx(1))*0.5,2))+ 2*(pow((u[0].dx(2) + u[2].dx(0))*0.5,2))+ 2*(pow((u[2].dx(1) + u[1].dx(2))*0.5,2))

        # if self.mesh.Dim == 2:
        #     TraceD2 =  pow(u[0].dx(0),2) + pow(u[1].dx(1),2) + 2*(pow((u[1].dx(0) + u[0].dx(1))*0.5,2)) 

        # gammadot = pow( abs(TraceD2) *0.5, 0.5)
        gammadot = self.gammaDot(u)
        

        iniPhi = max(project(localPhi*(phiInf-phi0)+phi0,V).vector().get_local())
        phi  = interpolate(Expression("iniPhi", degree=1,iniPhi = iniPhi), V)

        g_s = (phi/gammadot)*pow((abs(gammadot/phi-sigmay)/k),(1/nPow))/((phiInf-phi0)+(phi/gammadot)*pow((abs(gammadot/phi-sigmay)/k),(1/nPow)))
        c = 500
        H = 0.5*(1+tanh(c*(gammadot/phi-sigmay)))

        G = g_s - (phi-phi0)/(phiInf-phi0) 

        F = G * v * dx

        J = derivative(F,phi)

        bcs = []

        for sub in self.boundaries.inletBCs:
            bcs.append(
                DirichletBC(
                    V,
                    Constant(self.boundaries.Fluidityin*(phiInf-phi0)+phi0),
                    self.mesh.mf,
                    self.mesh.subdomains[sub],
                )
            )

        # Configurando o problema variacional
        problem = NonlinearVariationalProblem(F, phi,bcs, J)

        # Configurando o solucionador não linear
        solver = NonlinearVariationalSolver(problem)
        
        # Configurando os parâmetros do solucionador
        solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-06
        solver.parameters["newton_solver"]["relative_tolerance"] = 1e-07
        solver.parameters["newton_solver"]["maximum_iterations"] = 100
        solver.parameters['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        solver.solve()

        phi.rename("Phi Eq", "")
        Simulation_file = XDMFFile(f'phiEq.xdmf')
        Simulation_file.parameters["flush_output"] = True
        Simulation_file.parameters["functions_share_mesh"]= True
        Simulation_file.write(phi, 0.0)

        # ux = []
        # j = []
        # for i in np.linspace(-100e-6, 100e-6, 10):
        #     j.append(i)
        #     ux.append(phi(6e-3, i))
        # j

        return (phi-phi0)/(phiInf-phi0)

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


        a01 = (
            inner(
                self.TT(self.u, self.p, self.f),
                grad(self.v),
            )
        ) * self.mesh.dx()

        outletBCsIndex = tuple(
            self.mesh.subdomains[key]
            for key in self.boundaries.outletBCs
            if key in self.mesh.subdomains
        )
        # Outlet Pressure
        L01 = inner(dot(self.mesh.n , self.TT_direction(self.u, 0 , self.f,0)), self.v) * self.mesh.ds(outletBCsIndex)
        
        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(
                self.mesh.subdomains[key]
                for key in self.boundaries.inletBCs
                if key in self.mesh.subdomains
            )
            # Inlet Pressure
            L01 = L01+  inner(dot(self.mesh.n , self.TT_direction(self.u, 1 , self.f, 0)), self.v) * self.mesh.ds(inletBCsIndex)
            
        if self.boundaries.symmetryBCs!=None:
            symmetryBCsIndex = tuple(
                self.mesh.subdomains[key]
                for key in self.boundaries.symmetryBCs
                if key in self.mesh.subdomains
            )
            #Symmetry condition
            L01 = L01+  inner(dot(self.mesh.n , self.TT_direction(self.u, self.p , self.f, self.boundaries.symmetryAxis)), self.v) * self.mesh.ds(symmetryBCsIndex)
            
        # Mass Conservation(Continuity)
        a02 = (self.q * div(self.u)) * self.mesh.dx()
        L02 = 0

        # Fluidity
        if model=='newtonian' or model=='SMD':
            if model=='newtonian':
                phi = 1/self.fluid.k
            elif model=='SMD':
                phi = 1/self.etaSMD(self.fluid.k, self.fluid.nPow, self.u,1/self.fluid.phi0,1/self.fluid.phiInf)
        
            a03 = (
                (
                 self.normalized_fluidity(phi=phi,phi0=self.fluid.phi0,phiInf=self.fluid.phiInf)-self.f)/self.L
                * self.m
                * self.mesh.dx()
            )
            L03 = 0

        if model=='thixotropic':
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
                            *self.U
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
                            )
                        )
                    )
                ),
                (
                    -(
                        self.f
                        - dimensionless_phieq
                    )
                    / (
                        Tc
                       * self.U
                       )
                ),
            )


            a032 = (
                inner(
                    self.u,
                    grad(
                        self.f
                    ),
                )
            ) 

            a03 = a031 * self.m * self.mesh.dx()
            L03 = a032 * self.m * self.mesh.dx()

            #SUPG
            #Residual Strong form
            r3 = a031 - a032
            
            fnorm = sqrt(dot(self.f, self.f))
            delta = self.mesh.h/(2.0*fnorm)
            a03 +=delta*r3 *dot(self.u, grad(self.m))*self.mesh.dx()

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
        L01 = inner(dot(self.mesh.n , self.TT_direction(self.u, self.boundaries.Pout , (1 / (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0 )),0)), self.v) * self.mesh.ds(outletBCsIndex)
        
        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(
                self.mesh.subdomains[key]
                for key in self.boundaries.inletBCs
                if key in self.mesh.subdomains
            )
            # Inlet Pressure
            L01 = L01+  inner(dot(self.mesh.n , self.TT_direction(self.u, self.boundaries.Pin , (1 / (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0 )),0)), self.v) * self.mesh.ds(inletBCsIndex)

        if self.boundaries.symmetryBCs!=None:
            symmetryBCsIndex = tuple(
                self.mesh.subdomains[key]
                for key in self.boundaries.symmetryBCs
                if key in self.mesh.subdomains
            )
            #Symmetry condition
            L01 = L01+  inner(dot(self.mesh.n , self.TT_direction(self.u, self.p , (1 / (self.f *(self.fluid.phiInf - self.fluid.phi0) + self.fluid.phi0 )),self.boundaries.symmetryAxis)), self.v) * self.mesh.ds(symmetryBCsIndex)
            
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
        L03 = a033 * self.m * self.mesh.dx()

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