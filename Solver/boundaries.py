import sys
from dolfin import *
from mpi4py import MPI as pyMPI

sys.path.append("..")


class Boundaries:
    def __init__(self, mesh, UinVector_dim=None,UinMax_dim=None,Origin=None,R=None,Fluidityin = None ,Pout=0,inletBCs=["Inlet"],outletBCs=["Outlet"],noSlipBCs=["Wall"],symmetryBCs=None,symmetryAxis=None):
        ###########################################
        # inletCondition = 0 -> Constant inlet Pressure Condition
        # inletCondition = 1 -> Constant inlet X Velocity  Condition
        comm = MPI.comm_world
        
        self.mesh = mesh
        self.UinVector_dim = UinVector_dim
        self.UinMax_dim = UinMax_dim
        self.Origin = Origin
        self.R = R
        self.Fluidityin = Fluidityin
        self.inletBCs = inletBCs
        self.outletBCs = outletBCs
        self.noSlipBCs=noSlipBCs
        self.symmetryBCs=symmetryBCs
        
        # symmetryAxis
        # For x-axis 0
        # For y-axis 1
        # For z-axis 2
        self.symmetryAxis=symmetryAxis
        self.Pout = Pout
        self.bcs = []

        #! For future implementation: Inlet Pressure Condition
        # self.Pin = Pin
        # if Pin!=None:
        #     self.inletCondition=0
        if UinVector_dim != None or UinMax_dim!=None:
            self.inletCondition=1
        else:
            if comm.rank ==0:
                begin("\n***Velocity boundary conditions not provided.***\n***Please set either UinVector_dim or UinMax_dim to define inlet velocity.***")
            sys.exit()

        if self.mesh.Dim == 3:
            noSlipVector = Constant((0.0, 0.0, 0.0))
        elif self.mesh.Dim == 2:
            noSlipVector = Constant((0.0, 0.0))

        for sub in self.noSlipBCs:
            self.bcs.append(
                DirichletBC(
                    self.mesh.functionSpace.sub(0),
                    noSlipVector,
                    self.mesh.mf,
                    self.mesh.subdomains[sub],
                )
            )

        if self.symmetryBCs != None:
            for sub in self.symmetryBCs:
                self.bcs.append(
                    DirichletBC(
                        self.mesh.functionSpace.sub(0).sub(self.symmetryAxis),
                        Constant(0.0),
                        self.mesh.mf,
                        self.mesh.subdomains[sub],
                    )
                )

        if self.Fluidityin != None:
            for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(2),
                            Constant(self.Fluidityin),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )

        if self.inletCondition == 1:
            if UinVector_dim!=None:
                for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(0),
                            Constant(self.UinVector_dim),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )
            elif UinMax_dim!=None:
                if self.Origin==None or self.R==None:
                    if comm.rank ==0:
                        begin("\n***Origin or radius not provided.***\n***Please set Origin and R to define inlet velocity profile.***")
                    sys.exit()
                if self.mesh.Dim == 3:
                    u_D = Expression( ('UinMax_dim*(1-pow((pow( pow(x[0]-OriginX,2)+pow(x[1]-OriginY,2) ,0.5) )/R,2))','0','0') ,
                                     degree=2,UinMax_dim=UinMax_dim,R=R,OriginX=self.Origin[0],OriginY=self.Origin[1])
                elif self.mesh.Dim == 2:
                    u_D = Expression( ('UinMax_dim*(1-pow((pow( pow(x[0]-OriginX,2)+pow(x[1]-OriginY,2) ,0.5) )/R,2))','0') ,
                                     degree=2,UinMax_dim=UinMax_dim,R=R,OriginX=self.Origin[0],OriginY=self.Origin[1])
                for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(0),
                            u_D,
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )
        
    def change_parameter(self,Pin=None, UinVector=None,Fluidityin = None ,Pout=None,inletBCs=None,outletBCs=None,noSlipBCs=None):
        
        if Pin != None:
            self.Pin = Pin
        if UinVector!= None:
            self.UinVector = UinVector
        if Fluidityin!= None:
            self.Fluidityin = Fluidityin
        if inletBCs!= None:
            self.inletBCs = inletBCs
        if outletBCs!= None:
         self.outletBCs = outletBCs
        if noSlipBCs!= None:
            self.noSlipBCs=noSlipBCs
        if Pout!= None:
            self.Pout = Pout
        if Pin!=None:
            self.inletCondition=0
        elif UinVector != None:
            self.inletCondition=1


        if noSlipBCs!= None:
            self.bcs = []

            if self.mesh.Dim == 3:
                noSlipVector = Constant((0.0, 0.0, 0.0))
            elif self.mesh.Dim == 2:
                noSlipVector = Constant((0.0, 0.0))

            for sub in self.noSlipBCs:
                self.bcs.append(
                    DirichletBC(
                        self.mesh.functionSpace.sub(0),
                        noSlipVector,
                        self.mesh.mf,
                        self.mesh.subdomains[sub],
                    )
                )

        if self.Fluidityin != None:
            for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(2),
                            Constant(self.Fluidityin),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )

        if UinVector!= None:
            if self.inletCondition == 1:
                for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(0),
                            Constant(self.UinVector),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )