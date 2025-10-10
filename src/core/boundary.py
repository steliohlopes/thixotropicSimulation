import dolfin as df
from src.utils.config_parser import InletType

class Boundaries:
    def __init__(self, mesh, fluid, config):
        self.config = config.boundaries
        self.parameters = getattr(config, "parameters", None)

        V = mesh.function_space.sub(0)
        Q = mesh.function_space.sub(1)
        G = mesh.function_space.sub(2)
        
        self.bcs = []
        self.mesh = mesh
        self.fluid = fluid

    # --- compatibility attributes ---
    # TODO(legacy): remove this block when Problem/Solver consume new config directly
        self.inletBCs = tuple(self.config.inlet.subdomain)
        self.outletBCs = tuple(self.config.outlet)
        self.noSlipBCs = tuple(self.config.no_slip)
        self.Pout = self.parameters.outlet_pressure if self.parameters else 0.0
        self.R = self.config.inlet.radius
        origin = self.config.inlet.origin or (0.0,) * mesh.dimension
        if len(origin) < mesh.dimension:
            origin = origin + (0.0,) * (mesh.dimension - len(origin))
        self.Origin = origin
        self.UinVector_dim = None
        self.inletCondition = 0

        # Symmetry configuration for legacy solver
        if self.config.symmetry:
            axes = sorted(self.config.symmetry.items(), key=lambda kv: kv[0].value)
            symmetry_axes = [axis.value for axis, _ in axes]
            symmetry_sets = [tags for _, tags in axes]
            if len(symmetry_axes) == 1:
                self.symmetryAxis = symmetry_axes[0]
                self.symmetryBCs = symmetry_sets[0]
            else:
                self.symmetryAxis = symmetry_axes
                self.symmetryBCs = symmetry_sets
        else:
            self.symmetryAxis = None
            self.symmetryBCs = None
        
        # --- no slip ---
        zero_vec = df.Constant((0.0,)*mesh.dimension)
        
        for sub in self.config.no_slip:
            self.bcs.append(
                df.DirichletBC(
                    V,
                    zero_vec,
                    mesh.mf,
                    mesh.subdomains[sub],
                )
            )
        
        # --- symmetry ---
        for axis, tags in self.config.symmetry.items():
            for tag in tags:
                comp_space = V.sub(axis.value)
                self.bcs.append(
                    df.DirichletBC(
                        comp_space, 
                        df.Constant(0.0),
                        mesh.mf, mesh.subdomains[tag]
                    )
                )
                
        # --- inlet ---
        for subdomain in self.config.inlet.subdomain:
            if  self.config.inlet.type is InletType.VELOCITY_CONST:
                vec = df.Constant(self.config.inlet.value)
                self.UinVector_dim = self.config.inlet.value
                self.inletCondition = 1
                self.bcs.append(
                    df.DirichletBC(
                        V, 
                        vec,
                        mesh.mf, 
                        mesh.subdomains[subdomain]
                    )
                )
        
            elif self.config.inlet.type is InletType.FLOW_RATE:
                n = fluid.nPow if hasattr(fluid, "nPow") else (fluid.config.n or 1.0)
                axis = next(i for i,v in enumerate(self.config.inlet.value) if v != 0.0)
                flow_rate = self.config.inlet.value[axis]
                expr = self._powerlaw_inlet_expression(axis, flow_rate, self.R, n, self.Origin, mesh.dimension)
                self.inletCondition = 1
                self.bcs.append(
                    df.DirichletBC(
                        mesh.function_space.sub(0), 
                        expr,
                        mesh.mf, 
                        mesh.subdomains[subdomain]
                    )
                )
            elif self.config.inlet.type is InletType.PRESSURE:
                self.inletCondition = 0
                raise ValueError("Pressure inlets coming soon!  For now use 'velocity_const' or 'flow_rate'.")
            else:
                self.inletCondition = 0
                
        
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
            degree=degree, Q=Q, R=R, n=n, pi=df.pi,
            a0=origin[other[0]] if mesh_dim == 3 else origin[other],
            b0=origin[other[1]] if mesh_dim == 3 else 0.0,
        )
        return df.Expression(expr_tuple, **kwargs)
                                    
                