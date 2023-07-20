from ProblemInputs import Inputs
from PreProcessing.mesh import FiniteElementMesh


inputs = Inputs()
mesh = FiniteElementMesh(inputs)
mesh.msh2hdmf3D()
print(mesh)
