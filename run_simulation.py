from src.utils.config_parser import ConfigParser
from src.core.mesh import Mesh
from src.core.fluid import Fluid
from src.core.boundary import Boundaries 

json_file = 'inputs/base.json'

cfg = ConfigParser.from_file(json_file)

mesh = Mesh(cfg)
fluid = Fluid(cfg)
boundaries = Boundaries(mesh, fluid, cfg)
