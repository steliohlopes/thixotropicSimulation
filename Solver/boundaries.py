import sys
from dolfin import *

sys.path.append("..")


class Boundaries:
    def __init__(self, mesh, Pin=None, UinVector=None,Fluidityin = None ,Pout=0,inletBCs=["Inlet"],outletBCs=["Outlet"],noSlipBCs=["Wall"]):
        ###########################################
        # inletCondition = 0 -> Constant inlet Pressure Condition
        # inletCondition = 1 -> Constant inlet X Velocity  Condition
        # inletCondition = 2 -> fully developed flow inlet X Velocity  Condition

        self.mesh = mesh
        self.Pin = Pin
        self.UinVector = UinVector
        self.Fluidityin = Fluidityin
        self.inletBCs = inletBCs
        self.outletBCs = outletBCs
        self.noSlipBCs=noSlipBCs
        self.Pout = Pout
        self.bcs = []

        if Pin!=None:
            self.inletCondition=0
        elif UinVector != None:
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
            for sub in self.inletBCs:
                self.bcs.append(
                    DirichletBC(
                        self.mesh.functionSpace.sub(0),
                        Constant(self.UinVector),
                        self.mesh.mf,
                        self.mesh.subdomains[sub],
                    )
                )
