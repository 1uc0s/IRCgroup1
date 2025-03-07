#!/bin/bash
# Generated OpenFOAM setup commands for this HPC system

# Load required modules
module load mpi/intel-2018
module load openfoam/2.4.0

# Setup OpenFOAM environment
# Could not determine exact OpenFOAM setup path
# Try these approaches in your PBS file:
#
# Approach 1: Use path relative to project dir if available
if [ -n "$WM_PROJECT_DIR" ] && [ -f "$WM_PROJECT_DIR/etc/bashrc" ]; then
    source "$WM_PROJECT_DIR/etc/bashrc"
fi
#
# Approach 2: Check for common installation locations
for path in "/opt/openfoam240/etc/bashrc" "$HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc" "/usr/local/openfoam240/etc/bashrc"; do
    if [ -f "$path" ]; then
        source "$path"
        break
    fi
done

# Verify OpenFOAM environment
if command -v simpleFoam &> /dev/null; then
    echo "OpenFOAM environment set up successfully"
else
    echo "WARNING: OpenFOAM environment not properly configured"
fi
