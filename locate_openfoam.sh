#!/bin/bash
# locate_openfoam.sh - Script to identify OpenFOAM environment specifics on HPC
# This will help determine exactly how OpenFOAM is configured on your specific HPC

# Display colorful messages if terminal supports it
function echo_success() { echo "✓ $1"; }
function echo_info() { echo "ℹ $1"; }
function echo_warning() { echo "⚠ $1"; }
function echo_error() { echo "❌ $1"; }

echo_info "Locating OpenFOAM installation on HPC system..."

# Check module command
if command -v module &> /dev/null; then
    echo_success "Module command is available"
else
    echo_error "Module command not found!"
    exit 1
fi

# Check available modules
echo_info "Checking available OpenFOAM modules:"
module avail 2>&1 | grep -i foam

# Check Intel MPI modules
echo_info "Checking Intel MPI modules:"
module avail 2>&1 | grep -i "mpi/intel"

# Load OpenFOAM module
echo_info "Loading OpenFOAM 2.4.0 module..."
module load openfoam/2.4.0 2>&1

# Check environment variables after loading
echo_info "OpenFOAM environment variables:"
echo "FOAM_INST_DIR=$FOAM_INST_DIR"
echo "FOAM_APP=$FOAM_APP"
echo "FOAM_BASH=$FOAM_BASH"
echo "FOAM_LIBBIN=$FOAM_LIBBIN"

# Look for bashrc files
echo_info "Searching for OpenFOAM bashrc files..."
FOUND_FILES=$(find /opt /usr/local /apps $HOME -name "bashrc" -path "*OpenFOAM*2.4*" 2>/dev/null)

if [ -n "$FOUND_FILES" ]; then
    echo_success "Found OpenFOAM bashrc files:"
    echo "$FOUND_FILES"
else
    echo_warning "No OpenFOAM bashrc files found in standard locations."
fi

# Check for OpenFOAM executables
echo_info "Checking for OpenFOAM executables..."
if command -v simpleFoam &> /dev/null; then
    SIMPLEFOAM_PATH=$(which simpleFoam)
    echo_success "simpleFoam found at: $SIMPLEFOAM_PATH"
    echo_info "simpleFoam version info:"
    simpleFoam --help 2>&1 | head -3
else
    echo_error "simpleFoam not found in PATH!"
fi

if command -v blockMesh &> /dev/null; then
    BLOCKMESH_PATH=$(which blockMesh)
    echo_success "blockMesh found at: $BLOCKMESH_PATH"
else
    echo_error "blockMesh not found in PATH!"
fi

# Check OpenFOAM directory
if [ -n "$WM_PROJECT_DIR" ]; then
    echo_info "Checking OpenFOAM installation directory..."
    ls -la $WM_PROJECT_DIR/etc/
fi

# Check if OpenFOAM binaries are in PATH
echo_info "Checking PATH for OpenFOAM directories:"
echo $PATH | tr ':' '\n' | grep -i foam

# Output information for PBS file
echo_info "Creating OpenFOAM setup commands based on findings..."

OUTPUT_FILE="openfoam_setup_hpc.sh"
{
    echo "#!/bin/bash"
    echo "# Generated OpenFOAM setup commands for this HPC system"
    echo ""
    echo "# Load required modules"
    echo "module load mpi/intel-2018"
    echo "module load openfoam/2.4.0"
    echo ""
    echo "# Setup OpenFOAM environment"
    
    if [ -n "$FOAM_BASH" ]; then
        echo "# Using detected FOAM_BASH"
        echo "source $FOAM_BASH"
    elif [ -n "$FOUND_FILES" ]; then
        FIRST_FILE=$(echo "$FOUND_FILES" | head -1)
        echo "# Using detected bashrc file"
        echo "source $FIRST_FILE"
    else
        echo "# Could not determine exact OpenFOAM setup path"
        echo "# Try these approaches in your PBS file:"
        echo "#"
        echo "# Approach 1: Use path relative to project dir if available"
        echo "if [ -n \"\$WM_PROJECT_DIR\" ] && [ -f \"\$WM_PROJECT_DIR/etc/bashrc\" ]; then"
        echo "    source \"\$WM_PROJECT_DIR/etc/bashrc\""
        echo "fi"
        echo "#"
        echo "# Approach 2: Check for common installation locations"
        echo "for path in \"/opt/openfoam240/etc/bashrc\" \"\$HOME/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc\" \"/usr/local/openfoam240/etc/bashrc\"; do"
        echo "    if [ -f \"\$path\" ]; then"
        echo "        source \"\$path\""
        echo "        break"
        echo "    fi"
        echo "done"
    fi
    
    echo ""
    echo "# Verify OpenFOAM environment"
    echo "if command -v simpleFoam &> /dev/null; then"
    echo "    echo \"OpenFOAM environment set up successfully\""
    echo "else"
    echo "    echo \"WARNING: OpenFOAM environment not properly configured\""
    echo "fi"
} > $OUTPUT_FILE

chmod +x $OUTPUT_FILE
echo_success "Created setup script: $OUTPUT_FILE"
echo_info "Run this script to test OpenFOAM configuration or include its contents in your PBS file."

echo_info "Diagnosis complete. Please check the output for any issues."