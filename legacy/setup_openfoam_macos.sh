#!/bin/bash
# Setup script for OpenFOAM and Gmsh on macOS with mounted volumes

# Source OpenFOAM environment from mounted volume
if [ -d "/Volumes/OpenFOAM-v2412" ]; then
    # For OpenFOAM mounted as a volume
    export FOAM_INST_DIR="/Volumes"
    source "/Volumes/OpenFOAM-v2412/etc/bashrc"
    echo "✓ Sourced OpenFOAM from /Volumes/OpenFOAM-v2412"
else
    echo "❌ ERROR: OpenFOAM volume not found at /Volumes/OpenFOAM-v2412"
    echo "Please make sure the OpenFOAM volume is mounted and update the path if necessary."
    exit 1
fi

# Add Gmsh to PATH from mounted volume
if [ -d "/Volumes/gmsh-4.13.1-MacOSARM" ]; then
    # Check for the binary in typical locations within the volume
    if [ -f "/Volumes/gmsh-4.13.1-MacOSARM/bin/gmsh" ]; then
        export PATH="$PATH:/Volumes/gmsh-4.13.1-MacOSARM/bin"
        echo "✓ Added Gmsh to PATH from /Volumes/gmsh-4.13.1-MacOSARM/bin"
    elif [ -f "/Volumes/gmsh-4.13.1-MacOSARM/gmsh" ]; then
        export PATH="$PATH:/Volumes/gmsh-4.13.1-MacOSARM"
        echo "✓ Added Gmsh to PATH from /Volumes/gmsh-4.13.1-MacOSARM"
    else
        # Look for gmsh in the volume
        GMSH_PATH=$(find "/Volumes/gmsh-4.13.1-MacOSARM" -name "gmsh" -type f -perm +111 | head -1)
        if [ -n "$GMSH_PATH" ]; then
            export PATH="$PATH:$(dirname "$GMSH_PATH")"
            echo "✓ Added Gmsh to PATH from $(dirname "$GMSH_PATH")"
        else
            echo "❌ WARNING: Could not locate gmsh executable in the volume."
            echo "Please specify the full path to gmsh manually:"
            echo "export PATH=\$PATH:/path/to/gmsh/directory"
            exit 1
        fi
    fi
else
    echo "❌ ERROR: Gmsh volume not found at /Volumes/gmsh-4.13.1-MacOSARM"
    echo "Please make sure the Gmsh volume is mounted and update the path if necessary."
    exit 1
fi

# Verify that OpenFOAM and Gmsh are in PATH
echo "Checking OpenFOAM and Gmsh availability..."
if command -v simpleFoam >/dev/null 2>&1; then
    echo "✓ OpenFOAM found: $(which simpleFoam)"
else
    echo "❌ OpenFOAM commands not found in PATH. Environment may not be properly sourced."
    exit 1
fi

if command -v gmsh >/dev/null 2>&1; then
    echo "✓ Gmsh found: $(which gmsh)"
else
    echo "❌ Gmsh not found in PATH. Check if the executable is in the expected location."
    exit 1
fi

echo "✅ Environment setup completed successfully!"