#!/bin/bash
# setup_environment.sh - Unified script to set up OpenFOAM environment for aorta simulations

# Display colorful messages
function echo_success() { echo -e "\033[0;32m✓ $1\033[0m"; }
function echo_info() { echo -e "\033[0;34mℹ $1\033[0m"; }
function echo_warning() { echo -e "\033[0;33m⚠ $1\033[0m"; }
function echo_error() { echo -e "\033[0;31m❌ $1\033[0m"; }

echo_info "Setting up environment for OpenFOAM and Gmsh..."

# Detect operating system
OS="$(uname)"
echo_info "Detected operating system: $OS"

# macOS specific setup
if [ "$OS" == "Darwin" ]; then
    echo_info "Setting up for macOS environment"
    
    # Check if OpenFOAM is mounted as a volume (macOS specific)
    if [ -d "/Volumes/OpenFOAM-v2412" ]; then
        export FOAM_INST_DIR="/Volumes"
        source "/Volumes/OpenFOAM-v2412/etc/bashrc"
        echo_success "Sourced OpenFOAM from /Volumes/OpenFOAM-v2412"
    # Check for standard macOS installation paths
    elif [ -d "$HOME/OpenFOAM/OpenFOAM-v2412" ]; then
        export FOAM_INST_DIR="$HOME/OpenFOAM"
        source "$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc"
        echo_success "Sourced OpenFOAM from $HOME/OpenFOAM/OpenFOAM-v2412"
    else
        echo_error "OpenFOAM installation not found. Please make sure OpenFOAM is installed or mounted."
        echo_info "On macOS, OpenFOAM is typically mounted as a volume or installed in $HOME/OpenFOAM."
        exit 1
    fi

    # Find Gmsh - check multiple possible locations
    GMSH_PATHS=(
        "/Volumes/gmsh-4.13.1-MacOSARM/bin/gmsh"
        "/Volumes/gmsh-4.13.1-MacOSARM/gmsh"
        "/Volumes/gmsh-4.13.1-MacOSARM/Gmsh.app/Contents/MacOS/gmsh"
        "/Applications/Gmsh.app/Contents/MacOS/gmsh"
        "/opt/homebrew/bin/gmsh"
        "/usr/local/bin/gmsh"
    )

    for path in "${GMSH_PATHS[@]}"; do
        if [ -f "$path" ]; then
            GMSH_PATH=$(dirname "$path")
            export PATH="$PATH:$GMSH_PATH"
            echo_success "Added Gmsh to PATH from $GMSH_PATH"
            break
        fi
    done

    # If Gmsh is still not found, check mounted volumes more thoroughly
    if ! command -v gmsh >/dev/null 2>&1; then
        if [ -d "/Volumes/gmsh-4.13.1-MacOSARM" ]; then
            GMSH_PATH=$(find "/Volumes/gmsh-4.13.1-MacOSARM" -name "gmsh" -type f -perm +111 | head -1)
            if [ -n "$GMSH_PATH" ]; then
                export PATH="$PATH:$(dirname "$GMSH_PATH")"
                echo_success "Added Gmsh to PATH from $(dirname "$GMSH_PATH")"
            else
                echo_warning "Could not locate gmsh executable in the mounted volume."
            fi
        else
            echo_error "Gmsh not found. Please install Gmsh or mount the Gmsh volume."
            exit 1
        fi
    fi

# Linux setup
elif [ "$OS" == "Linux" ]; then
    echo_info "Setting up for Linux environment"
    
    # Check for OpenFOAM in standard Linux locations
    if [ -f "/opt/openfoam2412/etc/bashrc" ]; then
        source "/opt/openfoam2412/etc/bashrc"
        echo_success "Sourced OpenFOAM from /opt/openfoam2412"
    elif [ -f "$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc" ]; then
        source "$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc"
        echo_success "Sourced OpenFOAM from $HOME/OpenFOAM/OpenFOAM-v2412"
    else
        echo_error "OpenFOAM installation not found. Please make sure OpenFOAM is installed."
        exit 1
    fi
    
    # Check if Gmsh is installed
    if ! command -v gmsh >/dev/null 2>&1; then
        echo_error "Gmsh not found. Please install Gmsh (e.g., apt install gmsh or yum install gmsh)."
        exit 1
    else
        echo_success "Found Gmsh: $(which gmsh)"
    fi
else
    echo_error "Unsupported operating system: $OS"
    echo_info "This script currently supports macOS and Linux."
    exit 1
fi

# Verify that OpenFOAM and Gmsh are in PATH
echo_info "Checking OpenFOAM and Gmsh availability..."
if command -v simpleFoam >/dev/null 2>&1; then
    echo_success "OpenFOAM found: $(which simpleFoam)"
else
    echo_error "OpenFOAM commands not found in PATH. Environment may not be properly sourced."
    exit 1
fi

if command -v gmsh >/dev/null 2>&1; then
    echo_success "Gmsh found: $(which gmsh)"
else
    echo_error "Gmsh not found in PATH. Check if the executable is in the expected location."
    exit 1
fi

echo_success "Environment setup completed successfully!"