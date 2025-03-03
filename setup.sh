#!/bin/bash
# setup.sh - Comprehensive environment setup for aorta CFD simulation
# Creates Python virtual environment and verifies/installs dependencies

# Display colorful messages
function echo_success() { echo -e "\033[0;32m✓ $1\033[0m"; }
function echo_info() { echo -e "\033[0;34mℹ $1\033[0m"; }
function echo_warning() { echo -e "\033[0;33m⚠ $1\033[0m"; }
function echo_error() { echo -e "\033[0;31m❌ $1\033[0m"; }
function echo_section() { echo -e "\033[1;36m== $1 ==\033[0m"; }

# Ensure script exits on error
set -e

echo_section "Aorta CFD Simulation Environment Setup"
echo_info "Setting up a complete environment for OpenFOAM, Gmsh, and Python dependencies"

# Detect operating system
OS="$(uname)"
echo_info "Detected operating system: $OS"

# Create Python virtual environment
echo_section "Setting up Python virtual environment"

# Check for Python
if ! command -v python3 &>/dev/null; then
    echo_error "Python 3 is required but not installed."
    
    if [ "$OS" == "Darwin" ]; then
        echo_info "On macOS, you can install Python with: brew install python3"
    elif [ "$OS" == "Linux" ]; then
        echo_info "On Linux, you can install Python with: sudo apt-get install python3 python3-pip python3-venv"
    fi
    
    echo_info "Please install Python 3 and run this script again."
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo_info "Creating Python virtual environment in 'venv' directory..."
    python3 -m venv venv
    echo_success "Virtual environment created"
else
    echo_info "Using existing virtual environment in 'venv' directory"
fi

# Activate virtual environment
echo_info "Activating virtual environment..."
source venv/bin/activate

# Install required Python packages
echo_info "Installing required Python packages..."
pip install --upgrade pip
pip install numpy matplotlib

echo_success "Python environment configured successfully"

# Install OpenFOAM if needed
echo_section "Setting up OpenFOAM"

OPENFOAM_FOUND=false

# Check for OpenFOAM
if command -v simpleFoam &>/dev/null; then
    echo_success "OpenFOAM is already installed and available"
    OPENFOAM_FOUND=true
else
    echo_warning "OpenFOAM command 'simpleFoam' not found in PATH."
    
    # Try to source OpenFOAM from common locations
    if [ "$OS" == "Darwin" ]; then
        # Common macOS installation paths
        OPENFOAM_PATHS=(
            "/opt/openfoam2412/etc/bashrc"
            "/Volumes/OpenFOAM-v2412/etc/bashrc"
            "$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc"
        )
        
        for path in "${OPENFOAM_PATHS[@]}"; do
            if [ -f "$path" ]; then
                echo_info "Found OpenFOAM installation at: $path"
                echo_info "Sourcing OpenFOAM environment..."
                source "$path"
                OPENFOAM_FOUND=true
                break
            fi
        done
        
        if [ "$OPENFOAM_FOUND" = false ]; then
            echo_error "OpenFOAM not found on this system."
            echo_info "On macOS, OpenFOAM can be installed using:"
            echo_info "1. Download from https://openfoam.org/download/macos/"
            echo_info "2. Follow the installation instructions"
            echo_warning "Please install OpenFOAM v2412 and run this script again."
        fi
        
    elif [ "$OS" == "Linux" ]; then
        # Common Linux installation paths
        OPENFOAM_PATHS=(
            "/opt/openfoam2412/etc/bashrc"
            "/usr/lib/openfoam/openfoam2412/etc/bashrc"
            "$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc"
        )
        
        for path in "${OPENFOAM_PATHS[@]}"; do
            if [ -f "$path" ]; then
                echo_info "Found OpenFOAM installation at: $path"
                echo_info "Sourcing OpenFOAM environment..."
                source "$path"
                OPENFOAM_FOUND=true
                break
            fi
        done
        
        if [ "$OPENFOAM_FOUND" = false ]; then
            echo_error "OpenFOAM not found on this system."
            echo_info "On Ubuntu/Debian, OpenFOAM can be installed using:"
            echo_info "sudo add-apt-repository http://dl.openfoam.org/ubuntu"
            echo_info "sudo sh -c \"wget -O - https://dl.openfoam.org/gpg.key | apt-key add -\""
            echo_info "sudo apt-get update"
            echo_info "sudo apt-get install openfoam2412"
            echo_warning "Please install OpenFOAM v2412 and run this script again."
        fi
    fi
fi

# Add OpenFOAM sourcing to activation script for future sessions
if [ "$OPENFOAM_FOUND" = true ]; then
    echo_info "Adding OpenFOAM environment to virtual environment activation..."
    ACTIVATE_SCRIPT="venv/bin/activate"
    
    # Get the path to OpenFOAM bashrc
    FOAM_BASHRC=$(echo $FOAM_INST_DIR/OpenFOAM-*bashrc 2>/dev/null || echo "")
    
    if [ -n "$FOAM_BASHRC" ] && [ -f "$FOAM_BASHRC" ]; then
        # Check if OpenFOAM source is already in activate script
        if ! grep -q "OpenFOAM" "$ACTIVATE_SCRIPT"; then
            # Add OpenFOAM sourcing to the end of the activate script
            echo "# Source OpenFOAM environment" >> "$ACTIVATE_SCRIPT"
            echo "source $FOAM_BASHRC" >> "$ACTIVATE_SCRIPT"
            echo_success "OpenFOAM environment integrated with virtual environment"
        else
            echo_info "OpenFOAM environment already integrated with virtual environment"
        fi
    else
        echo_warning "Could not automatically integrate OpenFOAM with virtual environment"
        echo_info "You may need to source OpenFOAM manually in future sessions"
    fi
fi

# Install Gmsh if needed
echo_section "Setting up Gmsh"

GMSH_FOUND=false

# Check for Gmsh
if command -v gmsh &>/dev/null; then
    echo_success "Gmsh is already installed and available"
    GMSH_FOUND=true
else
    echo_warning "Gmsh command not found in PATH."
    
    if [ "$OS" == "Darwin" ]; then
        # Check common macOS Gmsh locations
        GMSH_PATHS=(
            "/Applications/Gmsh.app/Contents/MacOS/gmsh"
            "/opt/homebrew/bin/gmsh"
            "/usr/local/bin/gmsh"
            "/Volumes/gmsh-4.13.1-MacOSARM/bin/gmsh"
        )
        
        for path in "${GMSH_PATHS[@]}"; do
            if [ -f "$path" ]; then
                echo_info "Found Gmsh at: $path"
                echo_info "Adding Gmsh to PATH..."
                export PATH="$PATH:$(dirname "$path")"
                GMSH_FOUND=true
                
                # Add to activation script
                if ! grep -q "gmsh" "venv/bin/activate"; then
                    echo "# Add Gmsh to PATH" >> "venv/bin/activate"
                    echo "export PATH=\"\$PATH:$(dirname "$path")\"" >> "venv/bin/activate"
                    echo_success "Gmsh integrated with virtual environment"
                fi
                break
            fi
        done
        
        if [ "$GMSH_FOUND" = false ]; then
            echo_error "Gmsh not found on this system."
            echo_info "On macOS, Gmsh can be installed using:"
            echo_info "brew install gmsh"
            echo_info "Or downloaded from: https://gmsh.info/#Download"
            echo_warning "Please install Gmsh and run this script again."
        fi
        
    elif [ "$OS" == "Linux" ]; then
        echo_info "On Linux, attempting to install Gmsh automatically..."
        
        if command -v apt-get &>/dev/null; then
            echo_info "Detected apt package manager, installing Gmsh..."
            sudo apt-get update && sudo apt-get install -y gmsh
            GMSH_FOUND=true
        elif command -v yum &>/dev/null; then
            echo_info "Detected yum package manager, installing Gmsh..."
            sudo yum install -y gmsh
            GMSH_FOUND=true
        else
            echo_error "Could not automatically install Gmsh."
            echo_info "Please install Gmsh manually from: https://gmsh.info/#Download"
            echo_warning "After installation, run this script again."
        fi
    fi
fi

# Create a file to track environment status
if [ "$OPENFOAM_FOUND" = true ] && [ "$GMSH_FOUND" = true ]; then
    echo_section "Environment Verification"
    
    echo_info "Verifying OpenFOAM..."
    if command -v simpleFoam &>/dev/null; then
        echo_success "OpenFOAM verified: $(which simpleFoam)"
        echo "OpenFOAM_VERSION=$(simpleFoam --version 2>&1 | head -n 1)" > .env_status
    else
        echo_error "OpenFOAM verification failed"
    fi
    
    echo_info "Verifying Gmsh..."
    if command -v gmsh &>/dev/null; then
        echo_success "Gmsh verified: $(which gmsh)"
        echo "GMSH_VERSION=$(gmsh --version 2>&1 | head -n 1)" >> .env_status
    else
        echo_error "Gmsh verification failed"
    fi
    
    echo_info "Verifying Python environment..."
    echo "PYTHON_VERSION=$(python --version 2>&1)" >> .env_status
    echo "NUMPY_VERSION=$(python -c 'import numpy; print(f"numpy-{numpy.__version__}")' 2>&1)" >> .env_status
    echo "MATPLOTLIB_VERSION=$(python -c 'import matplotlib; print(f"matplotlib-{matplotlib.__version__}")' 2>&1)" >> .env_status
    
    echo_success "Environment setup completed successfully!"
    echo_info "To activate this environment in the future, run:"
    echo_info "source venv/bin/activate"
else
    echo_warning "Environment setup incomplete."
    echo_info "Please install the missing components and run this script again."
fi

# Copy setup_environment.sh contents for compatibility
if [ -f "setup_environment.sh" ]; then
    echo_info "Backing up original setup_environment.sh to setup_environment.sh.bak"
    cp setup_environment.sh setup_environment.sh.bak
fi

# Create an activation shortcut
echo_info "Creating a convenience activation script 'activate.sh'..."
cat > activate.sh << EOL
#!/bin/bash
# Quick activation script for the aorta CFD environment
source venv/bin/activate
echo -e "\033[0;32m✓ Aorta CFD environment activated\033[0m"
EOL
chmod +x activate.sh

echo_section "Next Steps"
echo_info "1. Activate the environment:   source venv/bin/activate"
echo_info "2. Run a simulation:           ./runStenosis.sh"
echo_info "3. Visualize results:          python plot_results.py"