import sys
from dolfin import *
from mpi4py import MPI as pyMPI

sys.path.append("..")


class Boundaries:
    def __init__(self, mesh,fluid=None, UinVector_dim=None,Origin=None,R=None,Pout=0,inletBCs=["Inlet"],outletBCs=["Outlet"],noSlipBCs=["Wall"],symmetryBCs=None,symmetryAxis=None):
        ###########################################
        # inletCondition = 0 -> Constant inlet Pressure Condition
        # inletCondition = 1 -> Constant inlet X Velocity  Condition
        comm = MPI.comm_world
        
        self.mesh = mesh
        self.fluid = fluid
        self.UinVector_dim = UinVector_dim
        self.Origin = Origin
        self.R = R
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

        #! For future implementation: Inlet Pressure Condition
        # self.Pin = Pin
        # if Pin!=None:
        #     self.inletCondition=0
        if self.UinVector_dim != None or self.fluid!=None:
            self.inletCondition=1
        else:
            if comm.rank ==0:
                begin("\n***Velocity boundary conditions not provided.***\n***Please set either UinVector_dim or fluid to define inlet velocity.***")
            sys.exit()

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
            if type(self.symmetryAxis) == int:
                for sub in self.symmetryBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(0).sub(self.symmetryAxis),
                            Constant(0.0),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )
            else:
                for i in range(len(self.symmetryAxis)):
                    for sub in self.symmetryBCs[i]:
                        self.bcs.append(
                            DirichletBC(
                                self.mesh.functionSpace.sub(0).sub(self.symmetryAxis[i]),
                                Constant(0.0),
                                self.mesh.mf,
                                self.mesh.subdomains[sub],
                            )
                        )

        if self.inletCondition == 1:
            if self.UinVector_dim!=None:
                for sub in self.inletBCs:
                    self.bcs.append(
                        DirichletBC(
                            self.mesh.functionSpace.sub(0),
                            Constant(self.UinVector_dim),
                            self.mesh.mf,
                            self.mesh.subdomains[sub],
                        )
                    )

            if self.Origin==None or self.R==None:
                if comm.rank ==0:
                    begin("\n***Origin or radius not provided.***\n***Please set Origin and R to define inlet velocity profile.***")
                sys.exit()
            if self.mesh.Dim == 3:
                '''
                Inlet Velocity power Law 
                u(y,z) = (Q/(pi*pow(R,3)))*((3*n+1)/(pow(R,(1/n))*(n+1))*(pow(R,((n+1)/n))-pow(pow(pow(y - y_0 ,2)+pow(z - z_0,2),0.5),((n+1)/n))))

                Inlet Velocity dimensionless power Law 
                u*(y,z) = 1-pow(pow(pow(y - y_0 ,2)+pow(z - z_0,2),0.5)/R,((n+1)/n))
                '''
                UinletDimVec = Expression( ('1-pow(pow(pow(x[1] - y_0 ,2)+pow(x[2] - z_0,2),0.5)/R,((n+1)/n))','0','0') ,
                                    degree=2,R=self.R,z_0=self.Origin[2],y_0=self.Origin[1],n=self.fluid.nPow)
            elif self.mesh.Dim == 2:
                # Dimensionless power-law inlet profile for 2D domains (mirrors 3D form)
                y_0 = self.Origin[1] if len(self.Origin) > 1 else 0.0
                UinletDimVec = Expression(
                    (
                        '1-pow(fabs(x[1] - y_0)/R,((n+1)/n))',
                        '0',
                    ),
                    degree=2,
                    R=self.R,
                    y_0=y_0,
                    n=self.fluid.nPow,
                )
            for sub in self.inletBCs:
                self.bcs.append(
                    DirichletBC(
                        self.mesh.functionSpace.sub(0),
                        UinletDimVec,
                        self.mesh.mf,
                        self.mesh.subdomains[sub],
                    )
                )
            
    def _powerlaw_inlet_expression(self, axis: int, Q: float, R: float, n: float, origin: tuple[float, ...], mesh_dim: int, degree: int = 2):
        """
        Return a vector Expression whose 'axis' component carries the
        volumetric-flow-rate profile and the other components are zero.
        """

        # ---------------- build scalar profile -----------------
        # choose the two coordinates transverse to the flow direction
        if mesh_dim == 3:
            other = [i for i in (0, 1, 2) if i != axis]
            coord_a, coord_b = other
            # radial distance r = sqrt( (x[coord_a]-a0)^2 + (x[coord_b]-b0)^2 )
            scalar = (
                "(Q/(pi*pow(R,3))) * ((3*n+1)/(pow(R,(1/n))*(n+1)))"
                "*( pow(R,((n+1)/n))"
                "  - pow(pow(pow(x[{a}] - a0, 2) + pow(x[{b}] - b0, 2),0.5), (n+1)/n) )"
            ).format(a=coord_a, b=coord_b)

        else:
            other = 1 if axis == 0 else 0
            scalar = (
                "(Q/(pi*pow(R,3))) * ((3*n+1)/(pow(R,(1/n))*(n+1)))"
                "*( pow(R,((n+1)/n))"
                "  - pow(pow(pow(x[{a}] - a0, 2),0.5), (n+1)/n) )"
            ).format(a=other)

        # ---------------- assemble vector ----------------------
        vec = ["0.0", "0.0", "0.0"]
        vec[axis] = scalar
        expr_tuple = tuple(vec[:mesh_dim])

        # ---------------- create Expression --------------------
        kwargs = dict(
            degree=degree, Q=Q, R=R, n=n, pi=pi,
            a0=origin[other[0]] if mesh_dim == 3 else origin[other],
            b0=origin[other[1]] if mesh_dim == 3 else 0.0,
        )
        return Expression(expr_tuple, **kwargs)