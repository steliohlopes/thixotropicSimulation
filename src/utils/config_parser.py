from typing import Any
from pathlib import Path
from dataclasses import dataclass, field
from enum import Enum
import json

class MeshDimension(Enum):
    D2 = 2
    D3 = 3
    
class FluidType(Enum):
    NEWTONIAN   = "newtonian"
    POWER_LAW   = "power_law"
    THIXOTROPIC = "thixotropic"
    SMD         = "smd"
    
class BoundaryType(Enum):
    INLET_VELOCITY      = "inlet_velocity"
    INLET_PRESSURE      = "inlet_pressure"
    OUTLET_PRESSURE     = "outlet_pressure"
    NO_SLIP             = "no_slip"
    SYMMETRY            = "symmetry"

class Axis(Enum):
    X = 0
    Y = 1
    Z = 2
    
class InletType(Enum):
    PRESSURE        = "pressure"
    VELOCITY_CONST  = "velocity_const"
    FLOW_RATE       = "flow_rate"           #vector, and in only one direction
    
        
@dataclass
class MeshConfig:
    file_path:          str
    dimension:          MeshDimension
    velocity_element:   tuple[str, int] | None = None
    pressure_element:   tuple[str, int] | None = None
    fluidity_element:   tuple[str, int] | None = None
    output_dir:         str = "results/mesh/"
    
@dataclass
class FluidConfig:
    type:       FluidType
    density:    float | None = None
    viscosity:  float | None = None
    k:          float | None = None
    n:          float | None = None
    phi0:       float | None = None
    phiInf:     float | None = None
    ta:         float | None = None
    tc:         float | None = None

    def __post_init__(self):
        if self.type is FluidType.NEWTONIAN:
            if self.viscosity is None:
                raise ValueError("Newtonian fluid needs 'viscosity'.")
        
        elif self.type is FluidType.POWER_LAW:
            if None in (self.k, self.n):
                raise ValueError("Power-law fluid needs both 'k' and 'n'.")
            
        elif self.type is FluidType.SMD:
            missing = [p for p in ("k","n","phi0", "phiInf", "ta", "tc")
                       if getattr(self, p) is None]
            if missing:
                raise ValueError(f"SMD fluid missing parameters: {', '.join(missing)}")

        else:
            raise ValueError(f"Unsupported FluidType {self.type}")
@dataclass
class InletBC:
    subdomain:  list
    type:       InletType
    value:      float | tuple[float, float] |tuple[float, float, float]
    origin:     tuple[float, float] | tuple[float, float, float] | None = None
    radius:     float | None = None
    
    def __post_init__(self):
        if isinstance(self.subdomain, str):
            self.subdomain = (self.subdomain,)
        if isinstance(self.value, list):
            self.value = tuple(self.value)
        if isinstance(self.origin, list):
            self.origin = tuple(self.origin)
        
        
            
@dataclass
class BoundaryConfig:
    inlet: list[InletBC]
    outlet: list[str]
    no_slip: list[str]
    symmetry: dict[Axis, list[str]] = field(default_factory=dict)
    
    def __post_init__(self):
        if isinstance(self.outlet, str):
            self.outlet = (self.outlet,)
        if isinstance(self.no_slip, str):
            self.no_slip = (self.no_slip,)
            
    
@dataclass
class SimulationParameters:
    model:              str = "newtonian"
    length_scale:       float = 1.0
    flow_rate:          float = 1.0
    outlet_pressure:    float = 0.0
    dimensional_output: bool = True
    save_checkpoint:    bool = False


@dataclass(slots=True)
class SimulationConfig:
    name:       str
    mesh:       MeshConfig
    fluid:      FluidConfig
    boundaries: BoundaryConfig
    parameters: SimulationParameters
    
class ConfigParser:
    
    @classmethod
    def from_file(cls, path: str | Path):
        with open(path, "r", encoding="utf-8") as f:
            raw_confing = json.load(f)
        return cls._build_config(raw_confing)
    
    @classmethod
    def _build_config(cls, raw: dict[str, Any]):
        
        return SimulationConfig(
            name        = raw.get("simulation_name", "simulation"),
            mesh        = cls._build_mesh(raw["mesh"]),
            fluid       = cls._build_fluid(raw["fluid"]),
            boundaries  = cls._build_boundaries(raw["boundaries"]),
            parameters  = cls._build_parameters(raw.get("simulation", {}))
            )
        
    @classmethod
    def _build_mesh(cls, mesh_data: dict[str, Any]):
        """Construct MeshConfig from raw data."""

        return MeshConfig(
            file_path           = mesh_data["file_path"],
            dimension           = MeshDimension(mesh_data["dimension"]),
            velocity_element    = mesh_data.get("velocity_element"),
            pressure_element    = mesh_data.get("pressure_element"),
            fluidity_element    = mesh_data.get("fluidity_element"),
            output_dir          = mesh_data.get("output_dir", "results/mesh/"),
        )

    @classmethod
    def _build_fluid(cls, fluid_data: dict[str, Any]):
        """Construct FluidConfig from raw data."""
        
        return FluidConfig(
            type        = FluidType(fluid_data["type"].lower()),
            density     = fluid_data.get("density"),
            viscosity   = fluid_data.get("viscosity"),
            k           = fluid_data.get("k"),
            n           = fluid_data.get("n"),
            phi0        = fluid_data.get("phi0"),
            phiInf      = fluid_data.get("phiInf"),
            ta          = fluid_data.get("ta"),
            tc          = fluid_data.get("tc")
        )
        
    @classmethod
    def _build_boundaries(cls, boundaries_data: dict[str, Any]):
        inlet = InletBC(
            subdomain   = boundaries_data.get("inlet").get("subdomain"),
            type        = InletType(boundaries_data.get("inlet").get("type")),
            value       = boundaries_data.get("inlet").get("value"),
            origin      = boundaries_data.get("inlet").get("origin"),
            radius      = boundaries_data.get("inlet").get("radius"),
        )
        
        sym_dict = { Axis[a.upper()]: v for a, v in boundaries_data.get("symmetry", {}).items() }
        
        return BoundaryConfig(
            inlet    = inlet,
            outlet   = boundaries_data.get("outlet"),
            no_slip  = boundaries_data.get("no_slip"),
            symmetry = sym_dict
        )
    
    @classmethod
    def _build_parameters(cls, simulation_data: dict[str, Any]):
        if not simulation_data:
            return SimulationParameters()

        return SimulationParameters(
            model              = simulation_data.get("model", "newtonian").lower(),
            length_scale       = simulation_data.get("length_scale", 1.0),
            flow_rate          = simulation_data.get("flow_rate", 1.0),
            outlet_pressure    = simulation_data.get("outlet_pressure", 0.0),
            dimensional_output = simulation_data.get("dimensional_output", True),
            save_checkpoint    = simulation_data.get("save_checkpoint", False),
        )
        
        
