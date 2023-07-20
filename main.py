from ProblemInputs import Inputs
from PreProcessing.mesh import Mesh


inputs = Inputs()
mesh = Mesh(inputs)
mesh.msh2hdmf3D()
print(mesh)
