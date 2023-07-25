import sys

sys.path.append("..")
from dolfin import *
from ProblemInputs import Inputs
from PreProcessing.mesh import FiniteElementMesh


class Boundaries:
    def __init__(self, inletCondition, input=Inputs(), mesh=FiniteElementMesh()):
        ###########################################
        # inletCondition = 0 -> Constant inlet Pressure Condition
        # inletCondition = 1 -> Constant inlet X Velocity  Condition
        # inletCondition = 2 -> fully developed flow inlet X Velocity  Condition

        # TODO Fazer para 2D verificando o tamanho do elemento
        self.inletCondition = inletCondition
        self.bcs = []

        ## No slip Boundaries
        for sub in input.noSlipBCs:
            self.bcs.append(
                DirichletBC(
                    mesh.functionSpace.sub(0),
                    Constant((0.0, 0.0, 0.0)),
                    mesh.mf,
                    mesh.subdomains[sub],
                )
            )

        if self.inletCondition == 0:
            input.pressureBC()
            self.Pin = input.Pin
            self.inletBCs = input.inletBCs
            # TODO Adicionar condições no Solver para Inlet Pressure

        elif self.inletCondition == 1:
            input.VelocityBC()
            self.Uin = input.Uin
            self.inletBCs = input.inletBCs
            for sub in self.inletBCs:
                self.bcs.append(
                    DirichletBC(
                        mesh.functionSpace.sub(0),
                        Constant((self.Uin, 0.0, 0.0)),
                        mesh.mf,
                        mesh.subdomains[sub],
                    )
                )
        
        

