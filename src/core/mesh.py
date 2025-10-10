from src.utils.config_parser import MeshDimension
from pathlib import Path
from typing import Optional
import meshio
import dolfin as df

class Mesh:
    def __init__(self, config):
        print("\n=== Mesh Initialization ===")
        
        self.config = config.mesh
        
        # Internal storage
        self.mesh_obj: Optional[df.Mesh] = None
        self.mf: Optional[df.MeshFunction] = None
        self.cf: Optional[df.MeshFunction] = None
        self.function_space: Optional[df.FunctionSpace] = None
        self.dx: Optional[df.Measure] = None
        self.ds: Optional[df.Measure] = None
        self.dimension: Optional[int] = None
        self.DoF: Optional[int] = None

        self._read_physical_names()

        Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)
        self._initialize_mesh()

    def _initialize_mesh(self):
        if not Path(self.config.file_path).exists():
            raise FileNotFoundError(f"Arquivo de malha não encontrado: {self.config.file_path}")
        
        mesh_path = Path(self.config.output_dir) / "mesh.xdmf"
        mf_path = Path(self.config.output_dir) / "mf.xdmf"
        
        if not (mesh_path.exists() and mf_path.exists()):
            self._process_mesh(mesh_path, mf_path)
                
        self._create_mesh_objects(mesh_path, mf_path)
        self._setup_function_space()
        
        print(f"\n✓ Mesh successfully initialized")
        print(f"  Cells: {self.mesh_obj.num_cells()}")
        print(f"  Vertices: {self.mesh_obj.num_vertices()}")
        print(f"  DoFs: {self.dof}")
    
    def _read_physical_names(self):
        physical_names = {}
        
        with open(self.config.file_path, "r") as fid:
            found = False
            finished = False
            
            for line in fid:
                if "$EndPhysicalNames" in line:
                    finished = True
                    break
                elif "$PhysicalNames" in line:
                    found = True
                elif found and not finished:
                    parts = line.split()
                    if len(parts) == 3:
                        name = parts[2].strip('"')
                        physical_names[name] = int(parts[1])
                        
        self.subdomains = physical_names
        
    def _process_mesh(self,mesh_path: Path, mf_path: Path):
        self.msh = meshio.read(self.config.file_path)
        
        # Determine cell types based on dimension
        if self.config.dimension == MeshDimension.D3:
            volume_cell_type = "tetra"
            face_cell_type = "triangle"
            dim_slice = slice(None)  # Use all coordinates (x,y,z)
        else:
            volume_cell_type = "triangle"
            face_cell_type = "line"
            dim_slice = slice(2)  # Use only x,y coordinates
        
        for key in self.msh.cell_data_dict["gmsh:physical"].keys():
            if key == face_cell_type:
                face_data =  self.msh.cell_data_dict["gmsh:physical"][key]
            elif key == volume_cell_type:
                volume_data = self.msh.cell_data_dict["gmsh:physical"][key]
        for cell in self.msh.cells:
            if cell.type == volume_cell_type:
                volume_cells = cell.data
            elif cell.type == face_cell_type:
                face_cells = cell.data
                
        # Create volume mesh
        volume_mesh = meshio.Mesh(
            points=self.msh.points[:, dim_slice],
            cells={volume_cell_type: volume_cells},
            cell_data={"physical": [volume_data]}
        )
        
        # Create face mesh
        face_mesh = meshio.Mesh(
            points=self.msh.points[:, dim_slice],
            cells=[(face_cell_type, face_cells)],
            cell_data={"physical": [face_data]}
        )
        
        # Write files
        meshio.write(mesh_path, volume_mesh)
        meshio.write(mf_path, face_mesh)

    def _create_mesh_objects(self,mesh_path: Path, mf_path: Path):
        self.mesh_obj = df.Mesh()
        with df.XDMFFile(str(mesh_path)) as infile:
            infile.read(self.mesh_obj)
        
        # Dimensão das facetas (dim-1)
        facet_dim = self.mesh_obj.topology().dim() - 1
        mvc = df.MeshValueCollection("size_t", self.mesh_obj, facet_dim)
        
        with df.XDMFFile(str(mf_path)) as infile:
            infile.read(mvc, "physical")
        self.mf = df.cpp.mesh.MeshFunctionSizet(self.mesh_obj, mvc)
        
        cell_dim = self.mesh_obj.topology().dim()
        mvc2 = df.MeshValueCollection("size_t", self.mesh_obj, cell_dim)
        
        with df.XDMFFile(str(mesh_path)) as infile:
            infile.read(mvc2, "physical")
        self.cf = df.cpp.mesh.MeshFunctionSizet(self.mesh_obj, mvc2)
        
        # Configura medidas e propriedades geométricas
        self._setup_mesh_properties()

    def _setup_mesh_properties(self):
        self.element_shape = self.mesh_obj.ufl_cell()
        self.dimension = self.mesh_obj.geometric_dimension()
        self.h = df.CellDiameter(self.mesh_obj)
        self.n = df.FacetNormal(self.mesh_obj)

        # Medidas para integração
        self.dx = df.Measure("dx", domain=self.mesh_obj, subdomain_data=self.cf)
        self.ds = df.Measure("ds", domain=self.mesh_obj, subdomain_data=self.mf)

    def _setup_function_space(self):
        elements = []
        if self.config.velocity_element: 
            vel_family, vel_order = self.config.velocity_element
            self.Uel = df.VectorElement(
                vel_family, self.element_shape, vel_order
            )
            elements.append(self.Uel)
        if self.config.pressure_element:
            pres_family, pres_order = self.config.pressure_element
            self.Pel = df.FiniteElement(
                pres_family, self.element_shape, pres_order
            )
            elements.append(self.Pel)
        if self.config.fluidity_element:
            fluid_family, fluid_order = self.config.fluidity_element
            self.Fel = df.FiniteElement(
                fluid_family, self.element_shape, fluid_order
            )
            elements.append(self.Fel)
            
        self.UPel = df.MixedElement(elements)
        self.function_space = df.FunctionSpace(self.mesh_obj, self.UPel)
        # TODO(legacy): remove this alias once solvers use `function_space` directly
        self.functionSpace = self.function_space
        self.dof = self.function_space.dim()
        # TODO(legacy): drop DoF alias after legacy solver refactor
        self.DoF = self.dof
        # TODO(legacy): replace meshObj usage with mesh_obj in downstream code
        self.meshObj = self.mesh_obj

        # Ensure attributes expected by legacy code exist
        if not hasattr(self, "Uel"):
            self.Uel = None
        if not hasattr(self, "Pel"):
            self.Pel = None
        if not hasattr(self, "Fel"):
            self.Fel = None

    # Legacy API compatibility -------------------------------------------------
    # TODO(legacy): remove once callers use `mesh.dx` property directly
    def dx(self, *args, **kwargs):
        return self._dx(*args, **kwargs) if args or kwargs else self._dx

    # TODO(legacy): remove once callers use `mesh.ds` property directly
    def ds(self, *args, **kwargs):
        return self._ds(*args, **kwargs) if args or kwargs else self._ds