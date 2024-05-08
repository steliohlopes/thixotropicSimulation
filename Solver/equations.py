from dolfin import *
import sys
import timeit
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI as pyMPI
import math
import csv
sys.path.append("..")

class Solver:    
    def __init__(self, problem,linearSolver = "mumps",nonlinearSolver = "newton",preconditioner = "default",absTol = 1e-9,relTol = 1e-10,maxIter = 30):

        self.problem = problem
        ## Solver Parameters
        self.linearSolver = linearSolver
        self.nonlinearSolver = nonlinearSolver
        self.absTol = absTol
        self.relTol = relTol
        self.maxIter = maxIter
        self.preconditioner = preconditioner
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
        prmU0['newton_solver']['krylov_solver']['absolute_tolerance'] = self.absTol
        prmU0['newton_solver']['krylov_solver']['relative_tolerance'] = self.relTol
        prmU0['newton_solver']['krylov_solver']['maximum_iterations'] = self.maxIter
        prmU0['newton_solver']['krylov_solver']['monitor_convergence'] = True
        prmU0['newton_solver']['krylov_solver']['report'] = True
        (no_iterations,converged) = self.solverU0.solve()
        #info(prmU,True)  #get full info on the parameters

        return self.problem.w
    
    def SaveSimulationData(self,filePath,fileName,dimensional=False):
        comm = MPI.comm_world
        num_procs = comm.Get_size()
        (self.u1, self.p1,self.f1) = self.problem.w.leaf_node().split()
        if dimensional:
            U = FunctionSpace(self.problem.mesh.meshObj, self.problem.mesh.Uel)
            self.u1 = project(self.u1*self.problem.U,U)
            Q = FunctionSpace(self.problem.mesh.meshObj, self.problem.mesh.Pel)
            self.p1 = project(self.p1*(self.problem.boundaries.Pin - self.problem.boundaries.Pout) + self.problem.boundaries.Pout,Q)
            V = FunctionSpace(self.problem.mesh.meshObj, self.problem.mesh.Fel)
            self.f1 = project(self.f1*(self.problem.fluid.phiInf - self.problem.fluid.phi0) + self.problem.fluid.phi0,V)
        self.u1.rename("Velocity Vector", "")
        self.p1.rename("Pressure", "")
        self.f1.rename("Fluidity", "")
        Simulation_file = XDMFFile(filePath+fileName+".xdmf")
        Simulation_file.parameters["flush_output"] = True
        Simulation_file.parameters["functions_share_mesh"]= True
        Simulation_file.write(self.u1, 0.0)
        Simulation_file.write(self.p1, 0.0)
        Simulation_file.write(self.f1, 0.0)
        Simulation_file.close()

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
    
    def velocity_plot(self,sweep_dict,velocity_coord,num_points,fileName):
        """
        Function to create a velocity plot by sweeping a line through a 2D or 3D simulation.

        Arguments:
            sweep_dict (dict): Dictionary specifying coordinates to sweep and their values.
                * Keys: coordinate representing the component (0-based indexing, e.g., 0 for "x", 1 for "y", 2 for "z").
                * Values: Lists for sweeping a coordinate or single values for fixed coordinates.
            velocity_coord (int): Index of the coordinate representing the velocity component (0-based indexing, e.g., 0 for "x", 1 for "y", 2 for "z").
            num_points (int): Number of points along the sweep line.
            fileName (str): Name of the output file (.csv)

        Returns:
            list: List of velocity values along the sweep line.
        """
        simulation_type = len(sweep_dict.keys())
        for key, value in sweep_dict.items():
            if isinstance(value, list):
                var_coord = key
                start_val = value[0]
                end_val = value[1]
                

        fix_values = sweep_dict.copy()
        del fix_values[var_coord]
        
        ux = []
        j = []

        for i in np.linspace(start_val*0.999, end_val*0.999, num_points):
            j.append(i)
          
            local_point_dict = fix_values.copy()
            local_point_dict[var_coord] = i
            if simulation_type ==3:
                local_point = [local_point_dict[0],local_point_dict[1],local_point_dict[2]]
            elif simulation_type ==2:
                local_point = [local_point_dict[0],local_point_dict[1]]
            ux.append(self.peval(self.u1,np.array(local_point))[velocity_coord])
        
        
        # plt.plot(ux,j,'r',)
        
        # plt.xlabel(r'$u$ [m/s]', fontsize=16)
        # plt.ylabel(r'r [m]', fontsize=16)
        # plt.title(f'Velocity Profile', fontsize=16)
        # plt.tight_layout()
        # plt.xlim(0, max(ux)*1.1) 
        # plt.ylim(start_val, end_val) 
        # plt.savefig(fileName)
        # plt.close()
        with open(fileName, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            for x, y in zip(ux, j):
                csv_writer.writerow([y, x])


    def viscoty_plot(self,arrayx,fileName):
        x = []
        eta = []
        phi0 = self.problem.fluid.phi0
        phiInf = self.problem.fluid.phiInf
        l = 4e-3
        r = 50e-6
        R = 100e-6

        for i in np.linspace(arrayx[0], arrayx[1], 200):
            x.append(i)
            if abs(i) <=  (l/2):
                f = r + (R-r) *math.sin(math.pi*abs(i)/l)
            else:
                f = R

            dimensionless_phi = self.peval(self.f1,np.array([i,f*0.99]))
            phi = dimensionless_phi*(phiInf-phi0)+phi0
            etaN = 1/phi
            eta.append(etaN)
        
        plt.plot(x,eta,'r',)
        
        plt.ylabel(r'$\eta$ [Pa.s]', fontsize=16)
        plt.xlabel(r'z [m]', fontsize=16)
        plt.xlim(arrayx[0],arrayx[1])
        plt.ylim(0,10)
        plt.title(f'Viscosity along wall', fontsize=16)
        plt.tight_layout()
        plt.savefig(fileName)
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