# CFD Simulation of Blood Flow in Atherosclerotic Aorta

This repository contains a self-contained framework for simulating blood flow through an aorta with atherosclerotic plaque using OpenFOAM and Gmsh.

## Overview

This simulation framework models blood flow through an aorta with varying degrees of stenosis (narrowing) to study hemodynamic changes associated with atherosclerosis. The simulation uses:

- **OpenFOAM v2412**: For computational fluid dynamics
- **Gmsh**: For parametric mesh generation
- **Python**: For visualization and post-processing
- **k-ω SST Turbulence Model**: For accurate flow modeling near walls

## Requirements

The setup script will check for and help install:

- **OpenFOAM v2412**
- **Gmsh 4.x or newer**
- **Python 3.x** with packages:
  - numpy
  - matplotlib
- **ParaView** (for visualization)

## Quick Start

### 1. One-Step Setup

Run the comprehensive setup script to create a virtual environment and verify all dependencies:

```bash
chmod +x setup.sh
./setup.sh
```

This script will:
1. Create a Python virtual environment
2. Install required Python packages
3. Check for OpenFOAM and Gmsh installations
4. Integrate all components into the virtual environment

### 2. Activate the Environment

After setup, activate the environment:

```bash
source venv/bin/activate
# Or use the convenience script
./activate.sh
```

### 3. Run a Simulation

To run a simulation with default settings (moderate stenosis, 50% occlusion):

```bash
./runStenosis.sh
```

To specify a different stenosis level:

```bash
./runStenosis.sh --level 0.3 --name mild
```

Available options:
- `--level`: Stenosis level from 0.0 (healthy) to 0.8 (severe)
- `--name`: A name for the simulation case
- `--mesh`: Mesh resolution (default: 0.5, smaller = finer)
- `--iterations`: Number of iterations (default: 2000)

### 4. Visualize Results

After simulation completes, visualize results using ParaView:

```bash
paraFoam -case aorta_simulation_moderate
```

Or generate plots using the Python script:

```bash
python plot_results.py
```

This will create plots in the `plots` directory, showing:
- Velocity profiles along centerline
- Pressure distribution
- Wall shear stress values

## Installation Details

### macOS

1. **OpenFOAM**:
   - Download from https://openfoam.org/download/macos/
   - Follow installation instructions for v2412

2. **Gmsh**:
   ```bash
   brew install gmsh
   ```
   
3. **Python & Dependencies**:
   ```bash
   brew install python3
   ```
   (Other dependencies installed by setup.sh)

### Linux (Ubuntu/Debian)

1. **OpenFOAM**:
   ```bash
   sudo add-apt-repository http://dl.openfoam.org/ubuntu
   sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
   sudo apt-get update
   sudo apt-get install openfoam2412
   ```

2. **Gmsh**:
   ```bash
   sudo apt-get install gmsh
   ```

3. **Python & Dependencies**:
   ```bash
   sudo apt-get install python3 python3-pip python3-venv
   ```
   (Other dependencies installed by setup.sh)

## File Structure

- `setup.sh` - Comprehensive environment setup script
- `activate.sh` - Convenience script to activate the environment
- `runStenosis.sh` - Script to run a simulation with configurable stenosis
- `plot_results.py` - Post-process and generate plots from results
- `aorta2d.geo` - Gmsh geometry template for the stenotic aorta
- `aortaCase/` - Template OpenFOAM case structure

## Physics Model

The simulation uses:

- **Fluid**: Blood (Newtonian approximation)
  - Density: 1060 kg/m³
  - Kinematic viscosity: 3.3×10⁻⁶ m²/s
- **Flow**: Steady-state incompressible (simpleFoam)
- **Turbulence**: k-ω SST model
- **Inlet**: Parabolic velocity profile (0.4 m/s average)
- **Outlet**: Zero pressure gradient
- **Walls**: No-slip condition

## Geometry Parameters

The aorta geometry is created parametrically with the following features:

- **Diameter**: 25mm (typical adult aorta)
- **Length**: 150mm segment
- **Stenosis**: Configurable constriction at the midpoint
  - Level: 0.0 (healthy) to 0.8 (severe)
  - Length: 20mm
  - Position: Center of vessel

## Extending the Model

This baseline model can be extended to:

1. **Pulsatile flow**: 
   - Replace simpleFoam with pimpleFoam
   - Define time-varying inlet velocity using codedFixedValue

2. **Patient-specific geometry**:
   - Import medical imaging data (CT/MRI)
   - Convert to STL and create mesh

3. **Non-Newtonian blood model**:
   - Modify transportProperties with a Carreau-Yasuda model

4. **Fluid-structure interaction**:
   - Include wall compliance
   - Use solidDisplacementFoam coupled with fluid solver

## Troubleshooting

If you encounter issues:

1. **Setup problems**:
   - Check `.env_status` file for version information
   - Verify OpenFOAM environment with `echo $FOAM_INST_DIR`

2. **Mesh generation problems**:
   - Adjust mesh size parameter with `--mesh 0.8` (larger = coarser)
   - Check for geometry errors in Gmsh

3. **Solver divergence**:
   - Reduce relaxation factors in fvSolution
   - Increase iterations with `--iterations 5000`

4. **Visualization issues**:
   - Make sure ParaView is compatible with your OpenFOAM version
   - Try running `reconstructPar` first if running in parallel

## License

This project is provided for educational and research purposes.

## Acknowledgments

- OpenFOAM: The Open Source CFD Toolbox (https://www.openfoam.com)
- Gmsh: A three-dimensional finite element mesh generator (https://gmsh.info)