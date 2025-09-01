from dolfin import *
import sys
from ufl_legacy import (real,conditional,tanh)
import timeit
import numpy as np

sys.path.append("..")


class Problem:
    def __init__(self, mesh, fluid, boundaries,L,Q):
        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries
        self.L = L
        self.Q = Q
        self.U = (self.Q/(pi*pow(self.boundaries.R,3)))*((3*self.fluid.nPow+1)/(pow(self.boundaries.R,(1/self.fluid.nPow))*(self.fluid.nPow+1))*(pow(self.boundaries.R,((self.fluid.nPow+1)/self.fluid.nPow))))

        self.start = timeit.default_timer()

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

    # # Deformation Tensor Dimensionless
    def DD(self, u_dim,direction=None):
        #Full Deformation Tensor
        # D = sym(as_tensor([ [        u[0].dx(0)           , (u[1].dx(0) + u[0].dx(1))*0.5 , (u[0].dx(2) + u[2].dx(0))*0.5 ],
        #                     [(u[1].dx(0) + u[0].dx(1))*0.5,         u[1].dx(1)            , (u[2].dx(1) + u[1].dx(2))*0.5 ],
        #                     [(u[0].dx(2) + u[2].dx(0))*0.5, (u[2].dx(1) + u[1].dx(2))*0.5 ,         u[2].dx(2)] ]))
        if direction==None:
            D = 0.5 * (nabla_grad(u_dim) + nabla_grad(u_dim).T)
        else:
            if self.mesh.Dim == 3:
                if direction==0:
                    D = sym(as_tensor([ [        0        ,     (u_dim[0].dx(1))*0.5         , (u_dim[0].dx(2))*0.5 ],
                                    [(u_dim[0].dx(1))*0.5,         u_dim[1].dx(1)            , (u_dim[2].dx(1) + u_dim[1].dx(2))*0.5 ],
                                    [(u_dim[0].dx(2))*0.5, (u_dim[2].dx(1) + u_dim[1].dx(2))*0.5 ,         u_dim[2].dx(2)] ]))
                    
                elif direction==1:
                    D = sym(as_tensor([ [        u_dim[0].dx(0)          , (u_dim[1].dx(0))*0.5 , (u_dim[0].dx(2) + u_dim[2].dx(0))*0.5 ],
                                [(u_dim[1].dx(0))*0.5,         0            , ( u_dim[1].dx(2))*0.5 ],
                                [(u_dim[0].dx(2) + u_dim[2].dx(0))*0.5, ( u_dim[1].dx(2))*0.5 ,         u_dim[2].dx(2)] ]))
                    
                elif direction==2:
                    D = sym(as_tensor([ [        u_dim[0].dx(0)           , (u_dim[1].dx(0) + u_dim[0].dx(1))*0.5 , ( u_dim[2].dx(0))*0.5 ],
                                    [(u_dim[1].dx(0) + u_dim[0].dx(1))*0.5,         u_dim[1].dx(1)            , (u_dim[2].dx(1))*0.5 ],
                                    [(      u_dim[2].dx(0))*0.5,            (u_dim[2].dx(1))*0.5          ,        0] ]))
                    
            elif self.mesh.Dim == 2:
                if direction==0:
                    D = sym(as_matrix([ [        0        ,     (u_dim[0].dx(1))*0.5      ],
                                    [(u_dim[0].dx(1))*0.5,         u_dim[1].dx(1)        ],]))
                    
                elif direction==1:
                    D = sym(as_matrix([ [u_dim[0].dx(0)      , (u_dim[1].dx(0))*0.5 ],
                                    [(u_dim[1].dx(0))*0.5,         0        ] ]))
                
        return D
    
    # Stress Tensor
    # def TT(self, u_dim, p_dim, phi_dim):
    #     phiInf = self.fluid.phiInf
    #     phi0 = self.fluid.phi0
    #     dphi = phiInf - phi0
        
    #     # Cartesian
    #     T = (2*self.L*dphi)/(phi_dim*dphi+phi0) * self.DD(u_dim) - (p_dim * Identity(len(u_dim)))
    #     return T
    
    def TT(self, u_dim, p_dim, phi_dim,direction=None):
        phiInf = self.fluid.phiInf
        phi0 = self.fluid.phi0
        dphi = phiInf - phi0
        # Cartesian
        T = (2*self.L*dphi)/(phi_dim*dphi+phi0) * self.DD(u_dim,direction) - (p_dim * Identity(len(u_dim)))
        
        return T

    def gammaDot(self, u_dim):
        """
        Calculates the dimensional magnitude of the strain rate tensor for a given velocity field.

        Args:
            u_dim : Dimensionless velocity field.

        Returns:
            float: The dimensional magnitude of the strain rate tensor.
        """
        DD = 0.5 * (nabla_grad(u_dim*self.U) + nabla_grad(u_dim*self.U).T)
        return pow(2 * inner(DD, DD), 0.5)
    
    def etaSMD(self, k, nPow, u_dim,eta0):
        """
        Calculates the viscosity using the Souza-Mendes-Dutra model.

        This function implements the viscosity model described in:

        Souza Mendes, Paulo R. and Dutra, Eduardo S. S..
        "Viscosity Function for Yield-Stress Liquids" Applied Rheology,
        vol. 14, no. 6, 2004, pp. 296-302. DOI:10.1515/arh-2004-0016

        Args:
            k (float): Consistency coefficient of the power-law model.
            nPow (float): Power-law exponent.
            u_dim: Dimensionless velocity field.
            eta0 (float): Reference viscosity at zero shear rate.

        Returns:
            viscosity field.
        """
        tauY=1e-5
        gammaDot = self.gammaDot(u_dim)
        return (1-exp(-eta0*gammaDot/tauY))*(tauY/gammaDot + k * pow(gammaDot, nPow - 1))
    
    def calculate_inlet_phieq(self,k,nPow,phi0, phiInf):
        V = FunctionSpace(self.mesh.meshObj, self.mesh.Fel)
        phi  = Function(V)
        v = TestFunction(V)
            
        Uinlet = Expression( ('(Q/(pi*pow(R,3)))*((3*n+1)/(pow(R,(1/n))*(n+1))*(pow(R,((n+1)/n))-pow(pow(pow(x[1] - y_0 ,2)+pow(x[2] - z_0,2),0.5),((n+1)/n))))') ,
                                    degree=2,R=self.boundaries.R,z_0=self.boundaries.Origin[2],y_0=self.boundaries.Origin[1],pi=pi,n=self.fluid.nPow,Q=self.Q)
        u = interpolate(Uinlet,V)
        if self.mesh.Dim == 3:
            TraceD2 = pow(u.dx(0),2) + pow((u.dx(1)),2) + pow(u.dx(2),2) 
        else:
            TraceD2 = pow(u.dx(0),2) + pow((u.dx(1)),2)
            
        gammadot = pow(abs(TraceD2) *0.5, 0.5)
        stress=pow(k*(gammadot),nPow)
        PhiIni = gammadot/(stress)
        phi=project(PhiIni, V)
        
        g_s = (phi/gammadot)*pow((abs(gammadot/phi)/k),(1/nPow))/((phiInf-phi0)+(phi/gammadot)*pow((abs(gammadot/phi)/k),(1/nPow)))
        c = 500
        H = 0.5*(1+tanh(c*(gammadot/phi)))

        G = g_s - (phi-phi0)/(phiInf-phi0) 

        F = G * v * dx

        J = derivative(F,phi)
        
        bcs = []
        
        # Configurando o problema variacional
        problem = NonlinearVariationalProblem(F, phi,bcs, J)

        # Configurando o solucionador não linear
        solver = NonlinearVariationalSolver(problem)
        
        # Configurando os parâmetros do solucionador
        solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10
        solver.parameters["newton_solver"]["relative_tolerance"] = 1e-10
        solver.parameters["newton_solver"]["maximum_iterations"] = 100
        solver.parameters['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        solver.solve()
        
        self.phiIn = interpolate(phi,V)
        self.phiIn_dim = interpolate(Expression('(phi-phi0)/(phiInf-phi0)',degree=2,phi=phi,phi0=phi0,phiInf=phiInf),V)
        
        return self.phiIn_dim

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
                    self.phiIn,
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

        return (phi-phi0)/(phiInf-phi0)

    def Tc(self):
        tc = 663
        return tc

    def Ta(self, dimensionless_phieq):
        ta = conditional(lt(dimensionless_phieq,1e-7),1e9,59.2 * (
            pow((1 - dimensionless_phieq), 1.1)
            / (pow(dimensionless_phieq, 0.4))
        ))

        return ta

    def S(self, dimensionless_phieq):
        s = conditional(lt(dimensionless_phieq,1e-3),1e9,(8 / (exp(dimensionless_phieq / 0.09) - 1)) + 1.2)
        return s

    
    def Equation(self,model ,wini=None,CheckPointFile=None):
        if CheckPointFile!=None:
            w_in =  XDMFFile(CheckPointFile)
            u_space = FunctionSpace(self.mesh.meshObj, self.mesh.Uel)
            p_space  = FunctionSpace(self.mesh.meshObj, self.mesh.Pel)
            f_space  = FunctionSpace(self.mesh.meshObj, self.mesh.Fel)

            u = Function(u_space)
            p = Function(p_space)
            f = Function(f_space)

            w_in.read_checkpoint(u, "u1")
            w_in.read_checkpoint(p, "p1")
            w_in.read_checkpoint(f, "f1")

            assign(self.w.sub(0), u)
            assign(self.w.sub(1), p)
            assign(self.w.sub(2), f)

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

        bcs = self.boundaries.bcs

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
        L01 = inner(dot(self.mesh.n , self.TT(self.u, 0 , self.f,0)), self.v) * self.mesh.ds(outletBCsIndex)
        
        #! For future implementation: Inlet Pressure Condition
        # if self.boundaries.inletCondition == 0:
        #     inletBCsIndex = tuple(
        #         self.mesh.subdomains[key]
        #         for key in self.boundaries.inletBCs
        #         if key in self.mesh.subdomains
        #     )
        #     # Inlet Pressure
        #     L01 = L01+  inner(dot(self.mesh.n , self.TT(self.u, 1 , self.f, 0)), self.v) * self.mesh.ds(inletBCsIndex)
            
        if self.boundaries.symmetryBCs!=None:
            if type(self.boundaries.symmetryAxis) == int:
                symmetryBCsIndex = tuple(
                    self.mesh.subdomains[key]
                    for key in self.boundaries.symmetryBCs
                    if key in self.mesh.subdomains
                )
                #Symmetry condition
                L01 = L01+  inner(dot(self.mesh.n , self.TT(self.u, self.p , self.f, self.boundaries.symmetryAxis)), self.v) * self.mesh.ds(symmetryBCsIndex)
            else:
                for i in range(len(self.boundaries.symmetryAxis)):
                    symmetryBCsIndex = tuple(
                        self.mesh.subdomains[key]
                        for key in self.boundaries.symmetryBCs[i]
                        if key in self.mesh.subdomains
                    )
                    #Symmetry condition
                    L01 = L01+  inner(dot(self.mesh.n , self.TT(self.u, self.p , self.f, self.boundaries.symmetryAxis[i])), self.v) * self.mesh.ds(symmetryBCsIndex)
            
        # Mass Conservation(Continuity)
        a02 = (self.q * div(self.u)) * self.mesh.dx()
        L02 = 0

        # Fluidity
        if model=='newtonian' or model=='SMD':
            if model=='newtonian':
                phi = 1/self.fluid.k
            elif model=='SMD':
                phi = 1/self.etaSMD(self.fluid.k, self.fluid.nPow, self.u,1/self.fluid.phi0)
        
            a03 = (
                (
                 ((phi - self.fluid.phi0) / (self.fluid.phiInf - self.fluid.phi0))-self.f)/self.L
                * self.m
                * self.mesh.dx()
            )
            L03 = 0

        if model=='thixotropic':
            phiIn_dim = self.calculate_inlet_phieq(k=self.fluid.k,nPow=self.fluid.nPow,phi0=self.fluid.phi0, phiInf=self.fluid.phiInf)

            for sub in self.boundaries.inletBCs:
                bcs.append(
                    DirichletBC(
                        self.mesh.functionSpace.sub(2),
                        phiIn_dim ,
                        self.mesh.mf,
                        self.mesh.subdomains[sub],
                    )
                )

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
            F0, self.w, bcs, J0
        )

        return self.problemU0