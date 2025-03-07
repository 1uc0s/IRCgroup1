#!/bin/bash
# Script to check OpenFOAM compatibility and alternatives on HPC

echo "==== OpenFOAM Compatibility Check ===="
echo "Checking CPU instruction set support..."

# Check if CPU info is available
if [ -f /proc/cpuinfo ]; then
    echo "CPU information found. Checking for required instructions..."
    
    # Check for each required instruction set
    for inst in mmx sse sse2 sse3 ssse3 sse4_1 sse4_2 avx popcnt; do
        if grep -q "$inst" /proc/cpuinfo; then
            echo "✓ $inst is supported"
        else
            echo "✗ $inst is NOT supported - may cause compatibility issues"
        fi
    done
else
    echo "Cannot access CPU information. Proceeding with caution."
fi

# Check available OpenFOAM modules
echo ""
echo "Checking available OpenFOAM versions..."
module avail openfoam 2>&1 | grep -i openfoam

echo ""
echo "==== Compatibility Test Plan ===="
echo "1. Try alternative OpenFOAM versions (recommended)"
echo "   - Look for older versions (e.g., 2.3.x) that may have fewer requirements"
echo "   - Check if there are any versions with 'compat' in the name"
echo ""
echo "2. Request specific nodes (if available)"
echo "   - Add to your PBS script: #PBS -l select=1:ncpus=4:mem=8gb:avx=true"
echo "   - Or consult your HPC documentation for available node types"
echo ""
echo "3. Use a pre-existing mesh"
echo "   - Convert a mesh created locally to OpenFOAM format"
echo "   - Then only run the solver on HPC (skipping mesh generation)"
echo ""
echo "4. Rebuild OpenFOAM with lower optimization"
echo "   - Contact your HPC administrator about this option"
echo ""

# Check if there are any other OpenFOAM versions available
echo "==== Finding Alternative Versions ===="
OTHER_VERSIONS=$(module avail 2>&1 | grep -i "foam" | grep -v "^--" | sort)

if [ -n "$OTHER_VERSIONS" ]; then
    echo "Found possible alternative versions:"
    echo "$OTHER_VERSIONS"
    
    # Look for older versions that might work
    OLDER_VERSION=$(echo "$OTHER_VERSIONS" | grep -i "openfoam/[12]" | head -1)
    if [ -n "$OLDER_VERSION" ]; then
        OLDER_VERSION=$(echo "$OLDER_VERSION" | awk '{print $1}')
        echo ""
        echo "Recommended first attempt: $OLDER_VERSION"
        echo "Add to your script: module load $OLDER_VERSION"
    fi
else
    echo "No alternative OpenFOAM versions found."
fi

# Create a simple test case to see if any version works
echo ""
echo "==== Creating a minimal test case ===="
TEST_DIR="openfoam_test"
mkdir -p $TEST_DIR/system
cd $TEST_DIR

# Create an extremely simple blockMeshDict
cat > system/blockMeshDict << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.0.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (0 0 0)
    (1 0 0)
    (1 1 0)
    (0 1 0)
    (0 0 1)
    (1 0 1)
    (1 1 1)
    (0 1 1)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (5 5 5) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    xMin
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    xMax
    {
        type patch;
        faces
        (
            (1 5 6 2)
        );
    }
    yMin
    {
        type patch;
        faces
        (
            (0 1 5 4)
        );
    }
    yMax
    {
        type patch;
        faces
        (
            (3 7 6 2)
        );
    }
    zMin
    {
        type patch;
        faces
        (
            (0 3 2 1)
        );
    }
    zMax
    {
        type patch;
        faces
        (
            (4 5 6 7)
        );
    }
);
EOL

# Generate a test script for testing multiple OpenFOAM versions
cat > test_versions.sh << EOL
#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -N openfoam_test
#PBS -j oe

cd \$PBS_O_WORKDIR/$TEST_DIR

echo "Testing all available OpenFOAM versions..."
module load tools/prod

# Get all OpenFOAM versions
VERSIONS=\$(module avail 2>&1 | grep -i "openfoam/" | grep -v "^--" | awk '{print \$1}')

if [ -z "\$VERSIONS" ]; then
    echo "No OpenFOAM versions found!"
    exit 1
fi

echo "Found OpenFOAM versions:"
echo "\$VERSIONS"
echo ""

for VERSION in \$VERSIONS; do
    echo "==================================================="
    echo "Testing \$VERSION"
    echo "==================================================="
    
    # Try to load the module
    module purge
    module load tools/prod
    module load \$VERSION
    
    if [ \$? -ne 0 ]; then
        echo "Failed to load \$VERSION, skipping..."
        continue
    fi
    
    # Try to run blockMesh
    echo "Running blockMesh with \$VERSION..."
    blockMesh > blockMesh_\$(echo \$VERSION | tr '/' '_').log 2>&1
    
    if [ \$? -eq 0 ]; then
        echo "SUCCESS: \$VERSION works on this system!"
        WORKING_VERSION=\$VERSION
        break
    else
        echo "FAILED: \$VERSION does not work on this system."
        echo "Error from log:"
        cat blockMesh_\$(echo \$VERSION | tr '/' '_').log | tail -5
    fi
done

if [ -n "\$WORKING_VERSION" ]; then
    echo ""
    echo "==================================================="
    echo "RECOMMENDATION: Use \$WORKING_VERSION in your scripts"
    echo "==================================================="
fi
EOL

chmod +x test_versions.sh

echo ""
echo "==== Next Steps ===="
echo "1. Run the test script to find a compatible OpenFOAM version:"
echo "   qsub test_versions.sh"
echo ""
echo "2. Look at the job output to see if any versions work."
echo ""
echo "3. If no versions work, consider creating the mesh locally and transferring it."
echo ""
echo "4. Contact HPC support for assistance with compatible software."

cd ..

echo ""
echo "Compatibility check complete!"