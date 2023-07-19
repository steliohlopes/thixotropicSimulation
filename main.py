from ProblemInputs import Inputs
from .PreProcessing.mesh import Mesh


inputs = Inputs()
mesh = Mesh(meshPath=inputs.meshPath, meshFile=inputs.meshFile)
