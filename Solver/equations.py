from dolfin import *
import sys
import timeit
sys.path.append("..")


class Solver:
    # Deformation Tensor
    def DD(self,u):
        # Cartesian
        D = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
        return D


    # Stress Tensor
    def TT(self,u, p, mu):
        # Cartesian
        T = 2 * mu * self.DD(u) - p * Identity(len(u))
        return T


    def gammaDot(self,u):
        return pow(2 * inner(self.DD(u), self.DD(u)), 0.5)


    def eta(self,k,nPow,u):
        eps=DOLFIN_EPS
        return k*pow(self.gammaDot(u)+eps,nPow-1)
    
    def __init__(self, mesh,fluid, boundaries):
        #TODO Criar classe Problem para receber o msh fluid e coundaries e definir a fisica do problema
        #TODO Retornar apenas o problemU0 para o SOlver apenas resolver o problema
        #TODO colocar as funçoes fora da classe, para dentro da classe


        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries

        #TODO Colocar como parametro do SOLVER, no parametro da classe
        ## Solver Parameters
        self.nonlinearSolver = "newton"
        self.absTol = 1e-9
        self.relTol = 1e-10
        self.maxIter = 30
        self.linearSolver = "mumps"

        ##### Functions
        ## Trial and Test function(s)
        self.dw = TrialFunction( self.mesh.functionSpace)
        (self.v, self.q) = TestFunctions( self.mesh.functionSpace) #Funçao peso
        self.w = Function( self.mesh.functionSpace)

        # Split into Velocity and Pressure
        if  self.mesh.Dim == 3:
            (self.u, self.p) = (as_vector((self.w[0], self.w[1], self.w[2])), self.w[3])
        elif  self.mesh.Dim == 2:
            (self.u, self.p) = (as_vector((self.w[0], self.w[1])), self.w[2])


    def NewtonianSolver(self,wini = None):
        self.start = timeit.default_timer()
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

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)
        #! Daqui pra cima mudar para Classe Problem
        solverU0 = NonlinearVariationalSolver(problemU0)
            # # Solver Parameters
        prmU0 = solverU0.parameters 
        prmU0['nonlinear_solver'] = self.nonlinearSolver
        prmU0['newton_solver']['absolute_tolerance'] = self.absTol
        prmU0['newton_solver']['relative_tolerance'] = self.relTol
        prmU0['newton_solver']['maximum_iterations'] = self.maxIter
        prmU0['newton_solver']['linear_solver'] = self.linearSolver
        prmU0['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        (no_iterations,converged) = solverU0.solve()
        self.stop = timeit.default_timer()

        return self.w
    
    def PowerLawSolver(self,wini = None):
        self.start = timeit.default_timer()
        if wini != None:
            self.w = wini
        
        a01 = (inner(self.TT(self.u,self.p,self.eta(self.fluid.k,self.fluid.nPow,self.u)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 
                                     
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(self.mesh.subdomains[self.boundaries.outletBCs[0]]) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(self.mesh.subdomains[self.boundaries.inletBCs[0]]) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)
        #! Aqui tambem
        solverU0 = NonlinearVariationalSolver(problemU0)
        # # Solver Parameters
        prmU0 = solverU0.parameters 
        prmU0['nonlinear_solver'] = self.nonlinearSolver
        prmU0['newton_solver']['absolute_tolerance'] = self.absTol
        prmU0['newton_solver']['relative_tolerance'] = self.relTol
        prmU0['newton_solver']['maximum_iterations'] = self.maxIter
        prmU0['newton_solver']['linear_solver'] = self.linearSolver
        prmU0['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        (no_iterations,converged) = solverU0.solve()
        self.stop = timeit.default_timer()

        return self.w
    
    def SaveSimulationData(self,fileName):
        (self.u1, self.p1) = self.w.leaf_node().split()
        self.u1.rename("Velocity Vector", "")
        self.p1.rename("Pressure", "")
        Simulation_file = XDMFFile(self.mesh.meshPath+fileName+".xdmf")
        Simulation_file.parameters["flush_output"] = True
        Simulation_file.parameters["functions_share_mesh"]= True
        Simulation_file.write(self.u1, 0.0)
        Simulation_file.write(self.p1, 0.0)

        total_time = self.stop - self.start
        # Output running time in a nice format.
        mins, secs = divmod(total_time, 60)
        hours, mins = divmod(mins, 60)

        log = open(self.mesh.meshPath+fileName+"_log.txt","w")
        log.write("rho="+str(self.fluid.rho)+
                  "\nk="+str(self.fluid.k)+
                  "\nnPow="+str(self.fluid.nPow)+
                  "\ntau0="+str(self.fluid.tau0)+
                  "\neta0="+str(self.fluid.eta0)+
                  "\netaInf="+str(self.fluid.etaInf)+
                  "\nts="+str(self.fluid.ts))
        log.write("\nTotal running time: %dh:%dmin:%ds \n" % (hours, mins, secs))
        #TODO Escrever numero de núcleos
        log.close()

        return