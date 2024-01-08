from dolfin import *
import sys
import timeit
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI as pyMPI
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
                  "\nTc="+str(self.problem.fluid.Tc)                  )
        log.write("\nTotal running time: %dh:%dmin:%ds" % (hours, mins, secs))
        log.write("\nNumber of processing cores utilized for the simulation: %d" % (num_procs))
        log.close()

        return
    
    def velocity_plot(self,R,xpoint,filePath):
        comm = MPI.comm_world
        # L length of pipe
        # R pipe radius
        nPow=self.problem.fluid.nPow
        k=self.problem.fluid.k
        Rloc = np.array([xpoint+DOLFIN_EPS_LARGE,0])
        Lloc = np.array([xpoint-DOLFIN_EPS_LARGE,0])
        dpdx= (self.peval(self.p1,Rloc)-self.peval(self.p1,Lloc))/(2*DOLFIN_EPS_LARGE)
        
        ux = []
        j = []
        uxPoiseuille = []
        uxPowerLaw=[]

        for i in np.linspace(-R*0.999, R*0.999, 200):
            j.append(i)
            ux.append(self.peval(self.u1,np.array([xpoint,i]))[0])
            uxPoiseuille.append(-(dpdx)*(1/(4*k))*(R**2-i**2))
            uxPowerLaw.append(nPow/(nPow+1)*(-dpdx/(2*k))**(1/nPow)*(R**((nPow+1)/nPow) -abs(i)**((nPow+1)/nPow)))

        
        plt.plot(
            uxPowerLaw/np.mean(uxPowerLaw), j,'g',
            uxPoiseuille/np.mean(uxPoiseuille), j,'b',
            ux/np.mean(ux),j,'r',)
        
        plt.xlabel(r'$\frac{u}{\bar{u}}$ [-]', fontsize=16)
        plt.ylabel(r'r [m]', fontsize=16)
        plt.title('Comparison of Velocity Profiles', fontsize=16)
        plt.legend([
            'PowerLaw Analytical',
            'Newtonian Analytical',
            'Thixotropic result'])
        plt.tight_layout()
        plt.savefig(filePath+'Result.png')
        plt.close()
    def mpi4py_comm(self,comm):
        '''Get mpi4py communicator'''
        try:
            return comm.tompi4py()
        except AttributeError:
            return comm

    
    def peval(self,f, x):
        '''Parallel synced eval'''
        try:
            yloc = f(x)
        except RuntimeError:
            yloc = np.inf*np.ones(f.value_shape())

        comm = self.mpi4py_comm(f.function_space().mesh().mpi_comm())
        yglob = np.zeros_like(yloc)
        comm.Allreduce(yloc, yglob, op=pyMPI.MIN)

        return yglob