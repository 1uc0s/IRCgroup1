# CFD Simulation of Blood Flow in Atherosclerotic Aorta

This repository contains files and scripts for simulating blood flow through an aorta with atherosclerotic plaque using OpenFOAM and Gmsh.

## Requirements

- OpenFOAM v2412
- Gmsh 4.x or newer
- Python 3.x with numpy and matplotlib
- ParaView (for visualization)

## File Structure

- `setup_openfoam_macos.sh` - Script to set up OpenFOAM environment on macOS
- `runSingleStenosis.sh` - Run a single case with configurable stenosis level
- `visualize_results.py` - Post-process and generate plots from results
- `.gitignore` - Ensures large meshes and output files aren't tracked by git

## Setup and Running

### 1. Environment Setup

First, ensure OpenFOAM and Gmsh are properly set up by running:

```bash
chmod +x setup_openfoam_macos.sh
source ./setup_openfoam_macos.sh
```

### 2. Run a Single Case

To run a simulation with default settings (moderate stenosis, 50% occlusion):

```bash
chmod +x runSingleStenosis.sh
./runSingleStenosis.sh
```

To specify a different stenosis level:

```bash
./runSingleStenosis.sh --level 0.3 --name mild
```

Available options:
- `--level` - Stenosis level from 0.0 (healthy) to 0.8 (severe), where the value represents the proportion of vessel diameter occluded
- `--name` - A name for the simulation case (for reference)

### 3. Visualize Results

After the simulation completes, you can visualize the results using ParaView:

```bash
paraFoam -case aorta_simulation
```

Or generate plots using the Python script:

```bash
python3 visualize_results.py
```

This will create various plots in the `plots` directory.

## Geometry and Mesh

The aorta geometry is created parametrically in Gmsh with the following features:

- Diameter: 25mm (typical adult aorta)
- Length: 150mm segment
- Stenosis: Parametrically controlled constriction at the midpoint
- 3D mesh: Created by extruding a 2D base geometry with a thin layer (quasi-2D)

## Physics Model

The simulation uses:

- Steady-state incompressible flow (simpleFoam)
- kOmegaSST turbulence model
- Blood properties:
  - Density: 1060 kg/m³
  - Kinematic viscosity: 3.3×10⁻⁶ m²/s
- Inlet velocity: 0.4 m/s (typical aortic blood flow)

## Extending the Model

This baseline model can be extended to:

1. **Pulsatile flow**: Replace simpleFoam with pimpleFoam and define time-varying inlet velocity
2. **Patient-specific geometry**: Import medical imaging data instead of parametric geometry
3. **Non-Newtonian blood model**: Modify transportProperties with a Carreau-Yasuda model
4. **Fluid-structure interaction**: Include wall compliance effects

## Troubleshooting

If you encounter issues:

1. **Mesh generation problems**: Adjust mesh size parameter in the .geo file
2. **Solver divergence**: Reduce relaxation factors in fvSolution
3. **Poor quality results**: Check y+ values to ensure proper mesh resolution near walls
4. **ParaView crashes**: Try reconstructPar first if running in parallel

## License

This project is provided for educational and research purposes.

## Acknowledgments

- OpenFOAM: The Open Source CFD Toolbox (https://www.openfoam.com)
- Gmsh: A three-dimensional finite element mesh generator (https://gmsh.info)