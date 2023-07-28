from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid


meshPath = "PreProcessing/Whistle/"
meshFile = "whistle"

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)
mesh.msh2hdmf3D()
mesh.createMeshObject3D

fluid = Fluid(
        rho=1000,
        k=0.1,
        nPow=0.572,
        tau0=6.21358,
        eta0=0.001,
        etaInf=64.1,
        ts=663
        )
