# Script to extract and plot data from OpenFOAM aorta simulations
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from pathlib import Path

# Base directory containing all simulation results
base_dir = "aorta_simulations"

# Stenosis cases to compare
cases = [
    {"name": "healthy", "dir": "stenosis_healthy", "color": "green"},
    {"name": "mild", "dir": "stenosis_mild", "color": "blue"},
    {"name": "moderate", "dir": "stenosis_moderate", "color": "orange"},
    {"name": "severe", "dir": "stenosis_severe", "color": "red"}
]

# Create output directory
os.makedirs("plots", exist_ok=True)

def extract_centerline_data(case_dir, field_name, time_dir="5000"):
    """Extract data along the centerline using sample utility"""
    
    # Create a sample dictionary file
    sample_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
        start   (0 0 0);
        end     (150 0 0);
        nPoints 500;
    }}
}}
"""
    
    # Write the sample dictionary
    with open(f"{case_dir}/system/sampleDict", "w") as f:
        f.write(sample_dict)
    
    # Run the sample utility
    os.system(f"postProcess -func sample -case {case_dir} -time {time_dir} > /dev/null 2>&1")
    
    # Read the sampled data
    data_file = Path(f"{case_dir}/{time_dir}/centerline_{field_name}_raw.xy")
    if not data_file.exists():
        print(f"Warning: No data file found at {data_file}")
        return None
    
    data = np.loadtxt(data_file)
    return data

# Plot velocity magnitude along centerline
def plot_velocity_centerline():
    plt.figure(figsize=(12, 6))
    
    for case in cases:
        case_dir = f"{base_dir}/{case['dir']}"
        
        # Get velocity data (U)
        u_data = extract_centerline_data(case_dir, "U")
        if u_data is None:
            continue
        
        # Calculate velocity magnitude
        x = u_data[:, 0]  # x-coordinate
        u_mag = np.sqrt(u_data[:, 1]**2 + u_data[:, 2]**2 + u_data[:, 3]**2)
        
        plt.plot(x, u_mag, label=case['name'].capitalize(), color=case['color'], linewidth=2)
    
    plt.xlabel("Position along aorta (mm)", fontsize=12)
    plt.ylabel("Velocity magnitude (m/s)", fontsize=12)
    plt.title("Blood Velocity along Aorta Centerline", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/velocity_centerline.png", dpi=300)
    plt.close()

# Plot pressure along centerline
def plot_pressure_centerline():
    plt.figure(figsize=(12, 6))
    
    for case in cases:
        case_dir = f"{base_dir}/{case['dir']}"
        
        # Get pressure data
        p_data = extract_centerline_data(case_dir, "p")
        if p_data is None:
            continue
        
        x = p_data[:, 0]  # x-coordinate
        p = p_data[:, 1]  # pressure
        
        plt.plot(x, p, label=case['name'].capitalize(), color=case['color'], linewidth=2)
    
    plt.xlabel("Position along aorta (mm)", fontsize=12)
    plt.ylabel("Pressure (m²/s²)", fontsize=12)
    plt.title("Pressure Distribution along Aorta Centerline", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/pressure_centerline.png", dpi=300)
    plt.close()

# Plot wall shear stress along the lower wall
def plot_wall_shear():
    plt.figure(figsize=(12, 6))
    
    for case in cases:
        case_dir = f"{base_dir}/{case['dir']}"
        
        # Create a sample dictionary file for the lower wall
        wall_sample_dict = f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
        
        # Write the sample dictionary
        with open(f"{case_dir}/system/wallSampleDict", "w") as f:
            f.write(wall_sample_dict)
        
        # Run the sample utility
        os.system(f"postProcess -dict system/wallSampleDict -case {case_dir} -time 5000 > /dev/null 2>&1")
        
        # Try to read the sampled data
        wss_file = None
        for file in Path(f"{case_dir}/5000").glob("lowerWall_*_raw.xy"):
            if "wallShearStress" in file.name:
                wss_file = file
                break
        
        if wss_file is None:
            print(f"Warning: No wall shear stress data found for {case['name']}")
            continue
        
        wss_data = np.loadtxt(wss_file)
        
        # Sort by x-coordinate
        wss_data = wss_data[wss_data[:, 0].argsort()]
        
        x = wss_data[:, 0]  # x-coordinate
        
        # Calculate WSS magnitude
        wss_mag = np.sqrt(wss_data[:, 1]**2 + wss_data[:, 2]**2 + wss_data[:, 3]**2)
        
        plt.plot(x, wss_mag, label=case['name'].capitalize(), color=case['color'], linewidth=2)
    
    plt.xlabel("Position along aorta (mm)", fontsize=12)
    plt.ylabel("Wall Shear Stress Magnitude (Pa)", fontsize=12)
    plt.title("Wall Shear Stress Distribution", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/wall_shear_stress.png", dpi=300)
    plt.close()

def main():
    print("Extracting and plotting results...")
    
    plot_velocity_centerline()
    plot_pressure_centerline()
    plot_wall_shear()
    
    print("Plotting completed. Results saved to 'plots' directory.")

if __name__ == "__main__":
    main()