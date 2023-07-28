from ProblemInputs import Inputs
from PreProcessing.mesh import FiniteElementMesh

meshPath = "PreProcessing/Whistle/"
meshFile = "whistle"

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)
mesh.msh2hdmf3D()
mesh.createMeshObject3D
