#!/usr/bin/env python3
# plot_results.py - Script to visualize results from aorta CFD simulations
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import argparse
from pathlib import Path
import subprocess
import re

def extract_centerline_data(case_dir, field_name, time_dir=None):
    """
    Extract data along the centerline of the aorta.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        field_name: Field to extract (e.g., 'U', 'p')
        time_dir: OpenFOAM time directory to use (default: latest time)
    
    Returns:
        numpy array with extracted data or None if extraction failed
    """
    case_path = Path(case_dir)
    
    # Find latest time directory if not specified
    if time_dir is None:
        time_dirs = [d for d in case_path.iterdir() if d.is_dir() and d.name.isdigit()]
        if not time_dirs:
            print(f"No time directories found in {case_dir}")
            return None
        time_dir = sorted(time_dirs, key=lambda x: float(x.name))[-1].name
    
    print(f"Using time directory: {time_dir}")
    
    # Create a sampling dictionary to extract data along the centerline
    sample_dict = f"""
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type            sets;
libs            (sampling);
interpolationScheme cellPoint;
setFormat       raw;

fields          ( {field_name} );

sets
{{
    centerline
    {{
        type    uniform;
        axis    distance;
        start   (0 0 0.5);
        end     (150 0 0.5);
        nPoints 500;
    }}
}}
"""
    
    # Make sure system directory exists
    system_dir = case_path / "system"
    system_dir.mkdir(exist_ok=True)
    
    # Write the sample dictionary
    sample_file = system_dir / "sampleDict"
    with open(sample_file, "w") as f:
        f.write(sample_dict)
    
    # Run the sample utility
    print(f"Sampling {field_name} along centerline...")
    subprocess.run(
        ["postProcess", "-func", "sample", "-case", str(case_path), "-time", time_dir],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Find the generated data file
    data_file = case_path / time_dir / f"centerline_{field_name}_raw.xy"
    
    if not data_file.exists():
        print(f"Warning: No data file found at {data_file}")
        return None
    
    # Read the data
    data = np.loadtxt(data_file)
    return data

def plot_velocity_centerline(case_dirs, labels=None, colors=None):
    """Plot velocity magnitude along the centerline for multiple cases"""
    if labels is None:
        labels = [os.path.basename(case_dir) for case_dir in case_dirs]
    
    if colors is None:
        cmap = get_cmap('viridis')
        colors = [cmap(i/len(case_dirs)) for i in range(len(case_dirs))]
    
    plt.figure(figsize=(12, 6))
    
    for i, case_dir in enumerate(case_dirs):
        # Get velocity data (U)
        u_data = extract_centerline_data(case_dir, "U")
        if u_data is None:
            print(f"Skipping {labels[i]} - no velocity data found")
            continue
        
        # Calculate velocity magnitude
        x = u_data[:, 0]  # x-coordinate
        u_mag = np.sqrt(u_data[:, 1]**2 + u_data[:, 2]**2 + u_data[:, 3]**2)
        
        plt.plot(x, u_mag, label=labels[i], color=colors[i], linewidth=2)
    
    plt.xlabel("Position along aorta (mm)", fontsize=12)
    plt.ylabel("Velocity magnitude (m/s)", fontsize=12)
    plt.title("Blood Velocity along Aorta Centerline", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    plt.savefig("plots/velocity_centerline.png", dpi=300)
    print("Saved plot to plots/velocity_centerline.png")
    
    return plt.gcf()

def plot_pressure_centerline(case_dirs, labels=None, colors=None):
    """Plot pressure along the centerline for multiple cases"""
    if labels is None:
        labels = [os.path.basename(case_dir) for case_dir in case_dirs]
    
    if colors is None:
        cmap = get_cmap('viridis')
        colors = [cmap(i/len(case_dirs)) for i in range(len(case_dirs))]
    
    plt.figure(figsize=(12, 6))
    
    for i, case_dir in enumerate(case_dirs):
        # Get pressure data
        p_data = extract_centerline_data(case_dir, "p")
        if p_data is None:
            print(f"Skipping {labels[i]} - no pressure data found")
            continue
        
        x = p_data[:, 0]  # x-coordinate
        p = p_data[:, 1]  # pressure
        
        plt.plot(x, p, label=labels[i], color=colors[i], linewidth=2)
    
    plt.xlabel("Position along aorta (mm)", fontsize=12)
    plt.ylabel("Pressure (m²/s²)", fontsize=12)
    plt.title("Pressure Distribution along Aorta Centerline", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    plt.savefig("plots/pressure_centerline.png", dpi=300)
    print("Saved plot to plots/pressure_centerline.png")
    
    return plt.gcf()

def extract_wall_shear_stress(case_dir, time_dir=None):
    """
    Extract wall shear stress data from the aorta wall.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        time_dir: OpenFOAM time directory to use (default: latest time)
    
    Returns:
        numpy array with extracted data or None if extraction failed
    """
    case_path = Path(case_dir)
    
    # Find latest time directory if not specified
    if time_dir is None:
        time_dirs = [d for d in case_path.iterdir() if d.is_dir() and d.name.isdigit()]
        if not time_dirs:
            print(f"No time directories found in {case_dir}")
            return None
        time_dir = sorted(time_dirs, key=lambda x: float(x.name))[-1].name
    
    # Create a sample dictionary for the lower wall
    wall_sample_dict = f"""
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type            sets;
libs            (sampling);
interpolationScheme cellPoint;
setFormat       raw;

fields          ( wallShearStress );

sets
{{
    lowerWall
    {{
        type    patchEdge;
        axis    x;
        patch   wall;
        pointDensity 300;
    }}
}}
"""
    # Make sure system directory exists
    system_dir = case_path / "system"
    system_dir.mkdir(exist_ok=True)
    
    # Write the sample dictionary
    with open(system_dir / "wallSampleDict", "w") as f:
        f.write(wall_sample_dict)
    
    # Run the sample utility
    print(f"Sampling wall shear stress...")
    subprocess.run(
        ["postProcess", "-dict", "system/wallSampleDict", "-case", str(case_path), "-time", time_dir],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Find the generated data file(s)
    data_file = None
    for file in (case_path / time_dir).glob("lowerWall_*_raw.xy"):
        if "wallShearStress" in file.name:
            data_file = file
            break
    
    if data_file is None:
        print(f"Warning: No wall shear stress data found for {case_dir}")
        return None
    
    # Read the data
    data = np.loadtxt(data_file)
    
    # Sort by x-coordinate
    data = data[data[:, 0].argsort()]
    
    return data

def plot_wall_shear_stress(case_dirs, labels=None, colors=None):
    """Plot wall shear stress for multiple cases"""
    if labels is None:
        labels = [os.path.basename(case_dir) for case_dir in case_dirs]
    
    if colors is None:
        cmap = get_cmap('viridis')
        colors = [cmap(i/len(case_dirs)) for i in range(len(case_dirs))]
    
    plt.figure(figsize=(12, 6))
    
    for i, case_dir in enumerate(case_dirs):
        # Get wall shear stress data
        wss_data = extract_wall_shear_stress(case_dir)
        if wss_data is None:
            print(f"Skipping {labels[i]} - no wall shear stress data found")
            continue
        
        x = wss_data[:, 0]  # x-coordinate
        
        # Calculate WSS magnitude
        wss_mag = np.sqrt(wss_data[:, 1]**2 + wss_data[:, 2]**2 + wss_data[:, 3]**2)
        
        plt.plot(x, wss_mag, label=labels[i], color=colors[i], linewidth=2)
    
    plt.xlabel("Position along aorta (mm)", fontsize=12)
    plt.ylabel("Wall Shear Stress Magnitude (Pa)", fontsize=12)
    plt.title("Wall Shear Stress Distribution", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    plt.savefig("plots/wall_shear_stress.png", dpi=300)
    print("Saved plot to plots/wall_shear_stress.png")
    
    return plt.gcf()

def find_simulation_dirs(base_dir=None, pattern=None):
    """Find simulation directories in the given base directory"""
    if base_dir is None:
        base_dir = os.getcwd()
    
    if pattern is None:
        pattern = r'aorta_simulation_.*'
    
    dirs = []
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path) and re.match(pattern, item):
            dirs.append(item_path)
    
    return sorted(dirs)

def extract_stenosis_level(case_dir):
    """Try to extract stenosis level from geo file in case directory"""
    geo_file = os.path.join(case_dir, "aorta.geo")
    if not os.path.exists(geo_file):
        return None
    
    with open(geo_file, 'r') as f:
        content = f.read()
    
    # Look for the stenosis level parameter
    match = re.search(r'stenosis_level\s*=\s*([\d.]+)', content)
    if match:
        return float(match.group(1))
    
    return None

def compare_stenosis_levels(case_dirs=None, labels=None, colors=None):
    """
    Create comparative plots for different stenosis levels.
    
    Args:
        case_dirs: List of case directories to compare
        labels: Labels for each case
        colors: Colors for each case
    """
    if case_dirs is None:
        case_dirs = find_simulation_dirs()
        
    if not case_dirs:
        print("No simulation directories found.")
        return
    
    # Extract stenosis levels if available
    if labels is None:
        labels = []
        for case_dir in case_dirs:
            level = extract_stenosis_level(case_dir)
            if level is not None:
                labels.append(f"Stenosis {level*100:.0f}%")
            else:
                labels.append(os.path.basename(case_dir))
    
    print(f"Found {len(case_dirs)} simulation directories:")
    for i, (case_dir, label) in enumerate(zip(case_dirs, labels)):
        print(f"{i+1}. {label}: {case_dir}")
    
    # Create plots directory
    os.makedirs("plots", exist_ok=True)
    
    # Plot velocity, pressure, and wall shear stress
    plot_velocity_centerline(case_dirs, labels, colors)
    plot_pressure_centerline(case_dirs, labels, colors)
    plot_wall_shear_stress(case_dirs, labels, colors)
    
    print("All plots have been saved to the 'plots' directory.")

def main():
    parser = argparse.ArgumentParser(description="Plot results from OpenFOAM aorta simulations")
    parser.add_argument("--dirs", nargs="+", help="Directories containing simulation results")
    parser.add_argument("--pattern", help="Pattern to match simulation directories")
    parser.add_argument("--labels", nargs="+", help="Labels for each simulation")
    parser.add_argument("--colors", nargs="+", help="Colors for each simulation")
    parser.add_argument("--velocity", action="store_true", help="Plot velocity along centerline")
    parser.add_argument("--pressure", action="store_true", help="Plot pressure along centerline")
    parser.add_argument("--wss", action="store_true", help="Plot wall shear stress")
    parser.add_argument("--all", action="store_true", help="Plot all data")
    
    args = parser.parse_args()
    
    # Find simulation directories
    if args.dirs:
        case_dirs = [os.path.abspath(d) for d in args.dirs]
    else:
        case_dirs = find_simulation_dirs(pattern=args.pattern)
    
    if not case_dirs:
        print("No simulation directories found.")
        return
    
    # Check if at least one plot type is requested
    if not (args.velocity or args.pressure or args.wss or args.all):
        args.all = True  # Default to plotting all
    
    # Plot requested data
    if args.all or args.velocity:
        plot_velocity_centerline(case_dirs, args.labels, args.colors)
    
    if args.all or args.pressure:
        plot_pressure_centerline(case_dirs, args.labels, args.colors)
    
    if args.all or args.wss:
        plot_wall_shear_stress(case_dirs, args.labels, args.colors)
    
    print("All plots have been saved to the 'plots' directory.")

if __name__ == "__main__":
    main()