import sys
sys.path.append("..")
from dolfin import *
from ProblemInputs import Inputs
from PreProcessing.mesh import FiniteElementMesh


class Boundaries:
    def __init__(self, inletCondition, input=Inputs(), mesh=FiniteElementMesh()):
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
