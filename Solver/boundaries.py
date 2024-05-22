import sys
from dolfin import *

sys.path.append("..")


class Boundaries:
    def __init__(self, mesh, Pin=None, UinVector=None,Umax_dim=None,Origin=None,R=None,Fluidityin = None ,Pout=0,inletBCs=["Inlet"],outletBCs=["Outlet"],noSlipBCs=["Wall"],symmetryBCs=None,symmetryAxis=None):
        ###########################################
        # inletCondition = 0 -> Constant inlet Pressure Condition
        # inletCondition = 1 -> Constant inlet X Velocity  Condition
        # inletCondition = 2 -> fully developed flow inlet X Velocity  Condition

        self.mesh = mesh
        self.Pin = Pin
        self.UinVector = UinVector
        self.Umax_dim = Umax_dim
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

        if Pin!=None:
            self.inletCondition=0
        elif UinVector != None or Umax_dim!=None:
            self.inletCondition=1

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
            if UinVector!=None:
                for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(0),
                            Constant(self.UinVector),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )
            elif Umax_dim!=None:
                if self.mesh.Dim == 3:
                    u_D = Expression( ('Umax_dim*(1-pow((pow( pow(x[0]-OriginX,2)+pow(x[1]-OriginY,2) ,0.5) )/R,2))','0','0') ,degree=2,Umax_dim=Umax_dim,R=R,OriginX=self.Origin[0],OriginY=self.Origin[1])
                elif self.mesh.Dim == 2:
                    u_D = Expression( ('Umax_dim*(1-pow((pow( pow(x[0]-OriginX,2)+pow(x[1]-OriginY,2) ,0.5) )/R,2))','0') ,degree=2,Umax_dim=Umax_dim,R=R,OriginX=self.Origin[0],OriginY=self.Origin[1])
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