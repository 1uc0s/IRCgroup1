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