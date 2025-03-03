#!/usr/bin/env python3
# plot_pulsatile_results.py - Enhanced script to visualize pulsatile flow results
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.cm import get_cmap
import argparse
from pathlib import Path
import subprocess
import re
import glob
from matplotlib.ticker import ScalarFormatter

def find_time_directories(case_dir):
    """Find all time directories in the OpenFOAM case directory"""
    time_dirs = []
    case_path = Path(case_dir)
    
    # Look for directories that are named as numbers
    for item in case_path.iterdir():
        if item.is_dir() and item.name.replace('.', '', 1).isdigit():
            time_dirs.append(item)
    
    # Sort by numerical value (not string)
    time_dirs.sort(key=lambda x: float(x.name))
    return time_dirs

def extract_data_for_time(case_dir, field_name, time_dir, sample_line="centerline"):
    """
    Extract field data for a specific time directory.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        field_name: Field to extract (e.g., 'U', 'p')
        time_dir: OpenFOAM time directory name
        sample_line: Name of sample line to extract data from
    
    Returns:
        numpy array with extracted data or None if extraction failed
    """
    case_path = Path(case_dir)
    
    # Build path to the sample data file
    sample_file = case_path / time_dir / f"{sample_line}_{field_name}_raw.xy"
    
    # If file doesn't exist, try to run the sample utility first
    if not sample_file.exists():
        print(f"Sample file for {field_name} not found for time {time_dir}. Running sample utility...")
        sample_dict_path = case_path / "system" / "sampleDict"
        
        # Check if sampleDict exists
        if not sample_dict_path.exists():
            print(f"No sampleDict found at {sample_dict_path}")
            return None
        
        # Run the sample utility for this time
        subprocess.run(
            ["postProcess", "-func", "sample", "-case", str(case_path), "-time", time_dir],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        # Check again if file was created
        if not sample_file.exists():
            print(f"Failed to create sample file: {sample_file}")
            return None
    
    # Read the data
    try:
        data = np.loadtxt(sample_file)
        return data
    except Exception as e:
        print(f"Error reading sample file {sample_file}: {e}")
        return None

def extract_velocity_magnitude(data):
    """Extract velocity magnitude from velocity vector data"""
    if data is None or data.shape[1] < 4:
        return None
    
    x = data[:, 0]  # x-coordinate
    u_mag = np.sqrt(data[:, 1]**2 + data[:, 2]**2 + data[:, 3]**2)
    return x, u_mag

def get_cycle_time(case_dir):
    """Try to extract cycle time from the velocity boundary condition"""
    case_path = Path(case_dir)
    u_file = case_path / "0" / "U"
    
    if not u_file.exists():
        return None
    
    # Look for cycleTime parameter in the U file
    with open(u_file, 'r') as f:
        content = f.read()
        
    match = re.search(r'cycleTime\s*=\s*([\d.]+)', content)
    if match:
        return float(match.group(1))
    
    # Alternative: look in controlDict for end time and try to infer
    control_dict = case_path / "system" / "controlDict"
    if control_dict.exists():
        with open(control_dict, 'r') as f:
            content = f.read()
        
        # Look for application type first
        app_match = re.search(r'application\s+(\w+)', content)
        if app_match and app_match.group(1) == "simpleFoam":
            return None  # Not a pulsatile simulation
            
        # Try to find end time and number of cycles
        end_match = re.search(r'endTime\s+(\d+\.?\d*)', content)
        if end_match:
            end_time = float(end_match.group(1))
            # Assume 5 cycles by default
            return end_time / 5
    
    return None

def plot_velocity_profile_animation(case_dir, output_file=None):
    """
    Create an animation of velocity profiles over time.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        output_file: Path to save the animation (mp4 format)
    """
    case_path = Path(case_dir)
    time_dirs = find_time_directories(case_dir)
    
    if not time_dirs:
        print(f"No time directories found in {case_dir}")
        return None
    
    print(f"Found {len(time_dirs)} time directories")
    
    # Get the cycle time if available
    cycle_time = get_cycle_time(case_dir)
    if cycle_time:
        print(f"Detected cardiac cycle time: {cycle_time} seconds")
    
    # Get initial data to setup plot
    initial_data = extract_data_for_time(case_dir, "U", time_dirs[0].name)
    if initial_data is None:
        print("Failed to extract initial velocity data")
        return None
    
    x, u_mag_initial = extract_velocity_magnitude(initial_data)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot(x, u_mag_initial, 'b-', linewidth=2)
    
    # Set up plot
    ax.set_xlabel("Position along aorta (mm)", fontsize=12)
    ax.set_ylabel("Velocity magnitude (m/s)", fontsize=12)
    ax.set_title("Blood Velocity along Aorta Centerline over Cardiac Cycle", fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Add a marker for stenosis location
    stenosis_pos = 75  # Default position (middle of the domain)
    ax.axvline(x=stenosis_pos, color='r', linestyle='--', alpha=0.5, label="Stenosis")
    
    # Add time text
    time_text = ax.text(0.02, 0.95, "Time: 0.000 s", transform=ax.transAxes)
    
    # Add legend
    ax.legend()
    
    # Get max velocity for consistent y-axis
    max_velocity = 0
    for time_dir in time_dirs[:min(10, len(time_dirs))]:  # Sample a few times
        data = extract_data_for_time(case_dir, "U", time_dir.name)
        if data is not None:
            _, u_mag = extract_velocity_magnitude(data)
            max_vel = np.max(u_mag)
            max_velocity = max(max_velocity, max_vel)
    
    # Add some margin
    max_velocity *= 1.1
    ax.set_ylim(0, max_velocity)
    
    # Update function for animation
    def update(frame):
        time_dir = time_dirs[frame]
        data = extract_data_for_time(case_dir, "U", time_dir.name)
        
        if data is not None:
            x, u_mag = extract_velocity_magnitude(data)
            line.set_data(x, u_mag)
            
            # Update time text
            time_value = float(time_dir.name)
            if cycle_time:
                cycle_phase = (time_value % cycle_time) / cycle_time
                time_text.set_text(f"Time: {time_value:.3f} s  (Phase: {cycle_phase:.2f})")
            else:
                time_text.set_text(f"Time: {time_value:.3f} s")
        
        return line, time_text
    
    # Create animation - use a subset of frames if there are too many
    num_frames = len(time_dirs)
    
    # If too many frames, sample them
    if num_frames > 100:
        step = num_frames // 100
        frames = range(0, num_frames, step)
    else:
        frames = range(num_frames)
    
    ani = FuncAnimation(fig, update, frames=frames, blit=True)
    
    if output_file:
        print(f"Saving animation to {output_file}")
        ani.save(output_file, writer='ffmpeg', fps=10, dpi=200)
    
    plt.tight_layout()
    return fig, ani

def plot_velocity_cycle_phases(case_dir, num_phases=8):
    """
    Plot velocity profiles at different phases of the cardiac cycle.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        num_phases: Number of phases to plot
    """
    case_path = Path(case_dir)
    time_dirs = find_time_directories(case_dir)
    
    if not time_dirs:
        print(f"No time directories found in {case_dir}")
        return None
    
    # Get the cycle time
    cycle_time = get_cycle_time(case_dir)
    if not cycle_time:
        print("Could not determine cardiac cycle time")
        return None
    
    print(f"Using cardiac cycle time: {cycle_time} seconds")
    
    # Group time directories by cycle phase
    phases = {}
    for time_dir in time_dirs:
        time_value = float(time_dir.name)
        phase = (time_value % cycle_time) / cycle_time
        phase_key = round(phase * num_phases) / num_phases
        
        if phase_key not in phases:
            phases[phase_key] = []
        
        phases[phase_key].append((time_dir, time_value))
    
    # Sort phases
    sorted_phases = sorted(phases.keys())
    
    # Select times closest to desired phases
    selected_times = []
    for phase in sorted_phases:
        if phases[phase]:
            # Use the last cycle for this phase
            times = sorted(phases[phase], key=lambda x: x[1])
            selected_times.append(times[-1])
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    cmap = get_cmap('viridis')
    colors = [cmap(i / len(selected_times)) for i in range(len(selected_times))]
    
    max_velocity = 0
    
    # Plot velocity profile for each phase
    for i, (time_dir, time_value) in enumerate(selected_times):
        data = extract_data_for_time(case_dir, "U", time_dir.name)
        if data is None:
            continue
        
        x, u_mag = extract_velocity_magnitude(data)
        phase = (time_value % cycle_time) / cycle_time
        
        ax.plot(x, u_mag, color=colors[i], linewidth=2, 
                label=f"Phase {phase:.2f} (t={time_value:.3f}s)")
        
        max_velocity = max(max_velocity, np.max(u_mag))
    
    # Add a marker for stenosis location
    stenosis_pos = 75  # Default position (middle of the domain)
    ax.axvline(x=stenosis_pos, color='r', linestyle='--', alpha=0.5)
    
    ax.set_xlabel("Position along aorta (mm)", fontsize=12)
    ax.set_ylabel("Velocity magnitude (m/s)", fontsize=12)
    ax.set_title("Blood Velocity Profiles at Different Phases of the Cardiac Cycle", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # Add some margin to y-axis
    ax.set_ylim(0, max_velocity * 1.1)
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    output_file = f"plots/velocity_phases_{os.path.basename(case_dir)}.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved phase plot to {output_file}")
    
    plt.tight_layout()
    return fig

def plot_wss_over_time(case_dir, num_locations=3):
    """
    Plot wall shear stress at specific locations over time.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        num_locations: Number of locations along the wall to track
    """
    case_path = Path(case_dir)
    time_dirs = find_time_directories(case_dir)
    
    if not time_dirs:
        print(f"No time directories found in {case_dir}")
        return None
    
    # Initialize data structures
    times = [float(time_dir.name) for time_dir in time_dirs]
    
    # First, extract geometry to determine interesting locations
    initial_data = extract_data_for_time(case_dir, "wallShearStress", time_dirs[0].name, "lowerWall")
    if initial_data is None:
        print("Failed to extract wall shear stress data")
        return None
    
    # Determine interesting locations (near stenosis)
    x_coords = initial_data[:, 0]
    x_range = max(x_coords) - min(x_coords)
    
    # Create positions array with focus on stenosis region
    center_pos = min(x_coords) + x_range / 2  # Stenosis center position
    stenosis_length = 20  # mm
    
    positions = [
        center_pos - stenosis_length,  # Before stenosis
        center_pos,                    # At stenosis
        center_pos + stenosis_length   # After stenosis
    ]
    
    # Find nearest indices for these positions
    position_indices = []
    for pos in positions:
        idx = np.argmin(np.abs(x_coords - pos))
        position_indices.append(idx)
    
    # Extract WSS over time at these positions
    wss_data = {pos: [] for pos in positions}
    
    for time_dir in time_dirs:
        data = extract_data_for_time(case_dir, "wallShearStress", time_dir.name, "lowerWall")
        if data is not None:
            for i, pos in enumerate(positions):
                idx = position_indices[i]
                wss_mag = np.sqrt(data[idx, 1]**2 + data[idx, 2]**2 + data[idx, 3]**2)
                wss_data[pos].append(wss_mag)
        else:
            # Fill with NaN if data not available
            for pos in positions:
                wss_data[pos].append(np.nan)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Get the cycle time if available
    cycle_time = get_cycle_time(case_dir)
    
    # If cycle time is available, add vertical lines for cycle boundaries
    if cycle_time:
        print(f"Detected cardiac cycle time: {cycle_time} seconds")
        
        max_time = max(times)
        cycle_boundaries = np.arange(0, max_time + cycle_time, cycle_time)
        
        for boundary in cycle_boundaries:
            if boundary <= max_time:
                ax.axvline(x=boundary, color='gray', linestyle='-', alpha=0.3)
    
    # Plot WSS for each location
    location_labels = ["Before Stenosis", "At Stenosis", "After Stenosis"]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, pos in enumerate(positions):
        ax.plot(times, wss_data[pos], label=f"{location_labels[i]} ({pos:.1f} mm)", 
                color=colors[i], linewidth=2)
    
    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel("Wall Shear Stress (Pa)", fontsize=12)
    ax.set_title("Wall Shear Stress at Different Locations Over Time", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # Use scientific notation for large numbers
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    output_file = f"plots/wss_over_time_{os.path.basename(case_dir)}.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved WSS plot to {output_file}")
    
    plt.tight_layout()
    return fig

def plot_pressure_drop_over_time(case_dir):
    """
    Plot pressure drop across the stenosis over time.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
    """
    case_path = Path(case_dir)
    time_dirs = find_time_directories(case_dir)
    
    if not time_dirs:
        print(f"No time directories found in {case_dir}")
        return None
    
    # Initialize data structures
    times = [float(time_dir.name) for time_dir in time_dirs]
    pressure_drops = []
    
    # First, extract geometry to determine stenosis location
    initial_data = extract_data_for_time(case_dir, "p", time_dirs[0].name)
    if initial_data is None:
        print("Failed to extract pressure data")
        return None
    
    # Determine stenosis location
    x_coords = initial_data[:, 0]
    x_range = max(x_coords) - min(x_coords)
    
    # Stenosis position
    center_pos = min(x_coords) + x_range / 2
    stenosis_length = 20  # mm
    
    # Positions for pressure measurements
    before_pos = center_pos - stenosis_length
    after_pos = center_pos + stenosis_length
    
    # Find indices for these positions
    before_idx = np.argmin(np.abs(x_coords - before_pos))
    after_idx = np.argmin(np.abs(x_coords - after_pos))
    
    # Extract pressure drop over time
    for time_dir in time_dirs:
        data = extract_data_for_time(case_dir, "p", time_dir.name)
        if data is not None:
            p_before = data[before_idx, 1]
            p_after = data[after_idx, 1]
            pressure_drop = p_before - p_after
            pressure_drops.append(pressure_drop)
        else:
            pressure_drops.append(np.nan)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Get the cycle time if available
    cycle_time = get_cycle_time(case_dir)
    
    # If cycle time is available, add vertical lines for cycle boundaries
    if cycle_time:
        print(f"Detected cardiac cycle time: {cycle_time} seconds")
        
        max_time = max(times)
        cycle_boundaries = np.arange(0, max_time + cycle_time, cycle_time)
        
        for boundary in cycle_boundaries:
            if boundary <= max_time:
                ax.axvline(x=boundary, color='gray', linestyle='-', alpha=0.3)
    
    # Plot pressure drop
    ax.plot(times, pressure_drops, linewidth=2, color='#d62728')
    
    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel("Pressure Drop (Pa)", fontsize=12)
    ax.set_title("Pressure Drop Across Stenosis Over Time", fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Use scientific notation for large numbers
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    output_file = f"plots/pressure_drop_{os.path.basename(case_dir)}.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved pressure drop plot to {output_file}")
    
    plt.tight_layout()
    return fig

def compare_velocity_profiles(case_dirs, labels=None, colors=None):
    """
    Compare velocity profiles at peak systole between different cases.
    
    Args:
        case_dirs: List of case directories to compare
        labels: Labels for each case
        colors: Colors for each case
    """
    if labels is None:
        labels = [os.path.basename(case_dir) for case_dir in case_dirs]
    
    if colors is None:
        cmap = get_cmap('viridis')
        colors = [cmap(i/len(case_dirs)) for i in range(len(case_dirs))]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    max_velocity = 0
    
    for i, case_dir in enumerate(case_dirs):
        # Find the time corresponding to peak systole
        time_dirs = find_time_directories(case_dir)
        
        if not time_dirs:
            print(f"No time directories found in {case_dir}")
            continue
        
        # Get the cycle time if available
        cycle_time = get_cycle_time(case_dir)
        if not cycle_time:
            print(f"Could not determine cardiac cycle time for {case_dir}")
            # Use latest time directory as fallback
            peak_time_dir = time_dirs[-1].name
        else:
            # Peak systole is typically around 0.2 of the cycle
            peak_phase = 0.2
            
            # Find the time directory closest to the desired phase
            # Take the last cycle for consistency
            max_time = float(time_dirs[-1].name)
            last_cycle_start = max(0, max_time - cycle_time)
            peak_time = last_cycle_start + peak_phase * cycle_time
            
            # Find closest time directory
            peak_time_dir = min(time_dirs, key=lambda x: abs(float(x.name) - peak_time)).name
        
        print(f"Using peak systole time {peak_time_dir} for {labels[i]}")
        
        # Extract velocity data at peak systole
        data = extract_data_for_time(case_dir, "U", peak_time_dir)
        if data is None:
            print(f"Failed to extract velocity data for {labels[i]}")
            continue
        
        x, u_mag = extract_velocity_magnitude(data)
        
        ax.plot(x, u_mag, label=labels[i], color=colors[i], linewidth=2)
        
        max_velocity = max(max_velocity, np.max(u_mag))
    
    # Add a marker for stenosis location
    stenosis_pos = 75  # Default position (middle of the domain)
    ax.axvline(x=stenosis_pos, color='r', linestyle='--', alpha=0.5, label="Stenosis")
    
    ax.set_xlabel("Position along aorta (mm)", fontsize=12)
    ax.set_ylabel("Velocity magnitude (m/s)", fontsize=12)
    ax.set_title("Peak Systole Velocity Profiles Comparison", fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # Add some margin to y-axis
    ax.set_ylim(0, max_velocity * 1.1)
    
    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)
    output_file = "plots/velocity_peak_comparison.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved velocity comparison plot to {output_file}")
    
    plt.tight_layout()
    return fig

def main():
    parser = argparse.ArgumentParser(description="Enhanced visualization for pulsatile OpenFOAM simulations")
    parser.add_argument("--case", help="Directory containing simulation results")
    parser.add_argument("--compare", nargs="+", help="Directories to compare")
    parser.add_argument("--labels", nargs="+", help="Labels for compared cases")
    parser.add_argument("--animate", action="store_true", help="Create animation of velocity profiles")
    parser.add_argument("--phases", action="store_true", help="Plot velocity at different cardiac phases")
    parser.add_argument("--wss", action="store_true", help="Plot wall shear stress over time")
    parser.add_argument("--pressure", action="store_true", help="Plot pressure drop over time")
    parser.add_argument("--all", action="store_true", help="Plot all available data")
    
    args = parser.parse_args()
    
    if args.compare:
        case_dirs = [os.path.abspath(d) for d in args.compare]
        compare_velocity_profiles(case_dirs, args.labels)
        return
    
    if args.case:
        case_dir = os.path.abspath(args.case)
    else:
        # Try to find the most recent simulation directory
        simulation_dirs = glob.glob("aorta_simulation_*")
        
        if not simulation_dirs:
            print("No simulation directories found. Please specify --case.")
            return
        
        # Use the most recently modified directory
        case_dir = max(simulation_dirs, key=os.path.getmtime)
        print(f"Using most recent simulation directory: {case_dir}")
    
    # Check if the case is a pulsatile simulation
    cycle_time = get_cycle_time(case_dir)
    
    if cycle_time is None:
        print(f"Case {case_dir} does not appear to be a pulsatile simulation.")
        print("For steady-state simulations, please use the regular plot_results.py script.")
        return
    
    # Determine which plots to generate
    if args.all or (not args.animate and not args.phases and not args.wss and not args.pressure):
        # Generate all plots
        if args.all or (len(find_time_directories(case_dir)) > 2):
            plot_velocity_cycle_phases(case_dir)
            plot_wss_over_time(case_dir)
            plot_pressure_drop_over_time(case_dir)
            
            # Animation is resource-intensive, only generate if explicitly requested
            if args.animate:
                output_file = f"plots/velocity_animation_{os.path.basename(case_dir)}.mp4"
                plot_velocity_profile_animation(case_dir, output_file)
    else:
        # Generate requested plots
        if args.animate:
            output_file = f"plots/velocity_animation_{os.path.basename(case_dir)}.mp4"
            plot_velocity_profile_animation(case_dir, output_file)
        
        if args.phases:
            plot_velocity_cycle_phases(case_dir)
        
        if args.wss:
            plot_wss_over_time(case_dir)
        
        if args.pressure:
            plot_pressure_drop_over_time(case_dir)
    
    print("All requested plots have been generated in the 'plots' directory.")

if __name__ == "__main__":
    main()