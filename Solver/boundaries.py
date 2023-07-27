import sys
from dolfin import *
from ProblemInputs import Inputs

sys.path.append("..")


class Boundaries:
    def __init__(self, inletCondition, mesh, input=Inputs()):
        ###########################################
        # inletCondition = 0 -> Constant inlet Pressure Condition
        # inletCondition = 1 -> Constant inlet X Velocity  Condition
        # inletCondition = 2 -> fully developed flow inlet X Velocity  Condition

        self.mesh = mesh
        self.inletCondition = inletCondition
        self.bcs = []

        ## No slip Boundaries
        for sub in input.noSlipBCs:
            if self.mesh.Dim == 3:
                noSlipVector = Constant((0.0, 0.0, 0.0))
            elif self.mesh.Dim == 2:
                noSlipVector = Constant((0.0, 0.0))

            self.bcs.append(
                DirichletBC(
                    self.mesh.functionSpace.sub(0),
                    noSlipVector,
                    self.mesh.mf,
                    self.mesh.subdomains[sub],
                )
            )

        if self.inletCondition == 0:
            input.pressureBC()
            self.Pin = input.Pin
            self.inletBCs = input.inletBCs

        elif self.inletCondition == 1:
            if self.mesh.Dim == 3:
                UinVector = Constant((input.Uin, 0.0, 0.0))
            elif self.mesh.Dim == 2:
                UinVector = Constant((input.Uin, 0.0))

            input.VelocityBC()
            self.Uin = input.Uin
            self.inletBCs = input.inletBCs
            for sub in self.inletBCs:
                self.bcs.append(
                    DirichletBC(
                        self.mesh.functionSpace.sub(0),
                        UinVector,
                        self.mesh.mf,
                        self.mesh.subdomains[sub],
                    )
                )
