from dolfin import *
import sys
import timeit
sys.path.append("..")

class Solver:    
    def __init__(self, problem,linearSolver = "mumps",nonlinearSolver = "newton",absTol = 1e-9,relTol = 1e-10,maxIter = 30):

        self.problem = problem
        ## Solver Parameters
        self.linearSolver = linearSolver
        self.nonlinearSolver = nonlinearSolver
        self.absTol = absTol
        self.relTol = relTol
        self.maxIter = maxIter
        self.start = timeit.default_timer()
    
    def SimulateEquation(self):
        
        self.solverU0 = NonlinearVariationalSolver(self.problem.problemU0)
        # Solver Parameters
        prmU0 = self.solverU0.parameters 
        prmU0['nonlinear_solver'] = self.nonlinearSolver
        prmU0['newton_solver']['absolute_tolerance'] = self.absTol
        prmU0['newton_solver']['relative_tolerance'] = self.relTol
        prmU0['newton_solver']['maximum_iterations'] = self.maxIter
        prmU0['newton_solver']['linear_solver'] = self.linearSolver
        prmU0['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        (no_iterations,converged) = self.solverU0.solve()

        return self.problem.w
    
    def SaveSimulationData(self,filePath,fileName):
        comm = MPI.comm_world
        num_procs = comm.Get_size()
        (self.u1, self.p1,self.f1) = self.problem.w.leaf_node().split()
        self.u1.rename("Velocity Vector", "")
        self.p1.rename("Pressure", "")
        self.f1.rename("Fluidity", "")
        Simulation_file = XDMFFile(filePath+fileName+".xdmf")
        Simulation_file.parameters["flush_output"] = True
        Simulation_file.parameters["functions_share_mesh"]= True
        Simulation_file.write(self.u1, 0.0)
        Simulation_file.write(self.p1, 0.0)
        Simulation_file.write(self.f1, 0.0)

        self.stop = timeit.default_timer()
        total_time = self.stop - self.start
        # Output running time in a nice format.
        mins, secs = divmod(total_time, 60)
        hours, mins = divmod(mins, 60)

        log = open(filePath+fileName+"_log.txt","w")
        log.write("rho="+str(self.problem.fluid.rho)+
                  "\nk="+str(self.problem.fluid.k)+
                  "\nnPow="+str(self.problem.fluid.nPow)+
                  "\ntau0="+str(self.problem.fluid.tau0)+
                  "\nphi0="+str(self.problem.fluid.phi0)+
                  "\nphiInf="+str(self.problem.fluid.phiInf)+
                  "\nTa="+str(self.problem.fluid.Ta)+
                  "\nphiTc="+str(self.problem.fluid.Tc)                  )
        log.write("\nTotal running time: %dh:%dmin:%ds" % (hours, mins, secs))
        log.write("\nNumber of processing cores utilized for the simulation: %d" % (num_procs))
        log.close()

        return