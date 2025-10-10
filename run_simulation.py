from pathlib import Path

from src.utils.config_parser import ConfigParser
from src.core.mesh import Mesh
from src.core.fluid import Fluid
from src.core.boundary import Boundaries
from src.core.problem import Problem
from src.core.equations import Solver


def main(config_path: str = "inputs/base.json") -> None:
	cfg = ConfigParser.from_file(config_path)

	mesh = Mesh(cfg)
	fluid = Fluid(cfg)
	boundaries = Boundaries(mesh, fluid, cfg)

	params = cfg.parameters
	model_map = {
		"smd": "SMD",
		"newtonian": "newtonian",
		"thixotropic": "thixotropic",
	}

	print("\n=== Problem setup ===")
	print(f"Model: {params.model}")
	print(f"Flow rate (Q): {params.flow_rate}")
	print(f"Length scale (L): {params.length_scale}")

	problem = Problem(
		mesh=mesh,
		fluid=fluid,
		boundaries=boundaries,
		L=params.length_scale,
		Q=params.flow_rate,
	)
	problem.Equation(model_map.get(params.model, params.model),CheckPointFile='results/PipeFlow2D/pipe_flow2D.h5')

	solver = Solver(problem)
	solver.SimulateEquation()

	results_dir = Path(mesh.config.output_dir)
	results_dir.mkdir(parents=True, exist_ok=True)
	file_path = str(results_dir)
	if not file_path.endswith("/"):
		file_path += "/"
	solver.SaveSimulationData(
		filePath=file_path,
		fileName=cfg.name,
		dimensional=params.dimensional_output,
		Checkpoint=params.save_checkpoint,
	)


if __name__ == "__main__":
	main()
