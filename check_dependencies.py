#!/usr/bin/env python3
# check_dependencies.py - Verify all dependencies are properly installed
import os
import sys
import subprocess
import platform
import shutil
from pathlib import Path

# ANSI color codes
class Colors:
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

def print_header(text):
    print(f"\n{Colors.BOLD}{Colors.BLUE}=== {text} ==={Colors.END}")

def print_success(text):
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")

def print_warning(text):
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")

def print_error(text):
    print(f"{Colors.RED}❌ {text}{Colors.END}")

def print_info(text):
    print(f"{Colors.BLUE}ℹ {text}{Colors.END}")

def run_command(command, error_msg=None):
    """Run a shell command and return output, or None if it fails"""
    try:
        result = subprocess.run(command, shell=True, check=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        if error_msg:
            print_error(error_msg)
        return None

def check_python():
    """Check Python version and packages"""
    print_header("Checking Python")
    
    # Check Python version
    py_version = platform.python_version()
    if py_version.startswith('3'):
        print_success(f"Python {py_version} found")
    else:
        print_error(f"Python 3.x required, found version {py_version}")
        return False
    
    # Check required packages
    required_packages = ['numpy', 'matplotlib']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            version = run_command(f"python -c 'import {package}; print({package}.__version__)'")
            print_success(f"{package} {version} installed")
        except ImportError:
            print_error(f"{package} not installed")
            missing_packages.append(package)
    
    if missing_packages:
        print_info("Missing packages can be installed with:")
        print_info(f"pip install {' '.join(missing_packages)}")
        return False
    
    return True

def check_openfoam():
    """Check OpenFOAM installation"""
    print_header("Checking OpenFOAM")
    
    # Check if OpenFOAM environment variables are set
    foam_dir = os.environ.get('FOAM_INST_DIR')
    
    if foam_dir:
        print_success(f"OpenFOAM environment found: {foam_dir}")
    else:
        print_warning("OpenFOAM environment not found ($FOAM_INST_DIR not set)")
        print_info("Checking for OpenFOAM executable...")
    
    # Check for simpleFoam executable
    simple_foam_path = shutil.which('simpleFoam')
    
    if simple_foam_path:
        # Get OpenFOAM version
        version_output = run_command("simpleFoam --version", "Failed to get OpenFOAM version")
        if version_output:
            version_line = version_output.split('\n')[0]
            print_success(f"simpleFoam found: {version_line}")
            
            # Check if version is 2412
            if '2412' not in version_output:
                print_warning(f"Expected OpenFOAM v2412, found: {version_line}")
                print_info("This simulation is designed for OpenFOAM v2412, other versions may cause issues")
        else:
            print_success("simpleFoam executable found (version unknown)")
        
        return True
    else:
        print_error("OpenFOAM (simpleFoam) not found in PATH")
        print_info("Make sure OpenFOAM is installed and sourced correctly")
        
        # Check common paths
        if platform.system() == "Darwin":  # macOS
            common_paths = [
                "/opt/openfoam2412/etc/bashrc",
                "/Volumes/OpenFOAM-v2412/etc/bashrc",
                os.path.expanduser("~/OpenFOAM/OpenFOAM-v2412/etc/bashrc")
            ]
        else:  # Linux
            common_paths = [
                "/opt/openfoam2412/etc/bashrc",
                "/usr/lib/openfoam/openfoam2412/etc/bashrc",
                os.path.expanduser("~/OpenFOAM/OpenFOAM-v2412/etc/bashrc")
            ]
        
        for path in common_paths:
            if os.path.exists(path):
                print_info(f"Found OpenFOAM installation at: {path}")
                print_info(f"Run: source {path}")
                break
        
        return False

def check_gmsh():
    """Check Gmsh installation"""
    print_header("Checking Gmsh")
    
    # Check for gmsh executable
    gmsh_path = shutil.which('gmsh')
    
    if gmsh_path:
        # Get Gmsh version
        version_output = run_command("gmsh --version", "Failed to get Gmsh version")
        if version_output:
            print_success(f"Gmsh found: {version_output}")
            
            # Parse version
            try:
                version = float(version_output.split(' ')[0])
                if version < 4.0:
                    print_warning(f"Gmsh version 4.0 or newer recommended, found: {version}")
            except:
                print_warning("Could not parse Gmsh version")
        else:
            print_success("Gmsh executable found (version unknown)")
        
        return True
    else:
        print_error("Gmsh not found in PATH")
        if platform.system() == "Darwin":  # macOS
            print_info("On macOS, install with: brew install gmsh")
        else:  # Linux
            print_info("On Ubuntu/Debian, install with: sudo apt-get install gmsh")
        
        return False

def check_paraview():
    """Check ParaView installation (optional)"""
    print_header("Checking ParaView (optional)")
    
    # Check for paraview executable
    paraview_path = shutil.which('paraview') or shutil.which('paraFoam')
    
    if paraview_path:
        if 'paraFoam' in paraview_path:
            print_success("paraFoam found")
        else:
            print_success("ParaView found")
        return True
    else:
        print_warning("ParaView not found (optional for visualization)")
        print_info("Download from: https://www.paraview.org/download/")
        return True  # Optional dependency

def check_virtual_env():
    """Check if running in a Python virtual environment"""
    print_header("Checking Virtual Environment")
    
    if os.environ.get('VIRTUAL_ENV'):
        venv_path = os.environ['VIRTUAL_ENV']
        print_success(f"Running in virtual environment: {venv_path}")
        return True
    else:
        print_warning("Not running in a virtual environment")
        print_info("For a self-contained setup, create and activate a virtual environment:")
        print_info("  python3 -m venv venv")
        print_info("  source venv/bin/activate (Linux/macOS)")
        print_info("  venv\\Scripts\\activate.bat (Windows)")
        return False

def check_simulation_files():
    """Check if all required simulation files exist"""
    print_header("Checking Simulation Files")
    
    required_files = [
        "aorta2d.geo",
        "runStenosis.sh",
        "plot_results.py",
        "setup.sh"
    ]
    
    required_dirs = [
        "aortaCase",
        "aortaCase/0",
        "aortaCase/system",
        "aortaCase/constant"
    ]
    
    missing_files = []
    missing_dirs = []
    
    for file in required_files:
        if os.path.isfile(file):
            print_success(f"Found file: {file}")
        else:
            print_error(f"Missing file: {file}")
            missing_files.append(file)
    
    for directory in required_dirs:
        if os.path.isdir(directory):
            print_success(f"Found directory: {directory}")
        else:
            print_error(f"Missing directory: {directory}")
            missing_dirs.append(directory)
    
    if missing_files or missing_dirs:
        print_warning("Some required files or directories are missing")
        return False
    
    return True

def write_status_file(status):
    """Write dependency status to a file"""
    with open('.dependency_status', 'w') as f:
        for component, is_ready in status.items():
            status_str = "READY" if is_ready else "MISSING"
            f.write(f"{component}={status_str}\n")
    
    print_info("Status written to .dependency_status")

def main():
    print_header("CFD Aorta Simulation Dependency Check")
    print_info("Checking all requirements for the CFD blood flow simulation...")
    
    status = {
        "python": check_python(),
        "openfoam": check_openfoam(),
        "gmsh": check_gmsh(),
        "paraview": check_paraview(),
        "virtual_env": check_virtual_env(),
        "simulation_files": check_simulation_files()
    }
    
    # Write status to file
    write_status_file(status)
    
    # Overall status
    print_header("Overall Status")
    
    ready = all([
        status["python"],
        status["openfoam"],
        status["gmsh"],
        status["simulation_files"]
    ])
    
    if ready:
        print_success("All essential dependencies are installed and ready!")
        print_info("You can now run the simulation with:")
        print_info("  ./runStenosis.sh")
    else:
        print_error("Some essential dependencies are missing")
        print_info("Please install the missing dependencies and run this script again")
    
    return 0 if ready else 1

if __name__ == "__main__":
    sys.exit(main())