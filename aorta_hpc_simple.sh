#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -N aorta_stenosis
#PBS -j oe

# =========================================================
# SIMPLIFIED SCRIPT FOR AORTA CFD ON HPC - NO CALCULATIONS
# =========================================================

# Go to the directory from which the job was submitted
cd $PBS_O_WORKDIR
echo "Working directory: $(pwd)"

# Load required modules
echo "Loading required modules..."
module load tools/prod
module load mpi/intel-2019
module load openfoam/2.4.0

# Source OpenFOAM environment
echo "Setting up OpenFOAM environment..."
if [ -f "/apps/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc" ]; then
    echo "Using OpenFOAM at /apps/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc"
    source "/apps/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc"
else
    echo "Could not find OpenFOAM environment file, but proceeding anyway..."
fi

# Verify OpenFOAM is available
if command -v simpleFoam &> /dev/null; then
    echo "OpenFOAM successfully loaded: $(which simpleFoam)"
    simpleFoam --version 2>&1 | head -1
else
    echo "ERROR: simpleFoam command not found. OpenFOAM not properly configured."
    exit 1
fi

# Set hardcoded values instead of calculations
STENOSIS_LEVEL="0.5"
CASE_NAME="moderate"
MESH_SIZE="0.5"
MAX_ITERATIONS="1000"

echo "======================================================="
echo "Simulation parameters:"
echo "  Stenosis level: $STENOSIS_LEVEL"
echo "  Case name: $CASE_NAME"
echo "  Mesh size: $MESH_SIZE"
echo "  Max iterations: $MAX_ITERATIONS"
echo "======================================================="

# Create case directory
case_dir="aorta_simulation_${CASE_NAME}"
echo "Creating simulation in: $case_dir"

# Create directories
if [ -d "$case_dir" ]; then
    echo "Directory $case_dir already exists. Cleaning up..."
    rm -rf "$case_dir"
fi

mkdir -p "$case_dir"
mkdir -p "$case_dir/0"
mkdir -p "$case_dir/constant"
mkdir -p "$case_dir/system"

# Create all OpenFOAM configuration files
echo "Creating OpenFOAM configuration files..."

# Create the controlDict
control_dict_file="$case_dir/system/controlDict"
cat > "$control_dict_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     simpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          1;

writeControl    timeStep;

writeInterval   200;

purgeWrite      3;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   8;

runTimeModifiable true;

functions
{
    wallShearStress
    {
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        writeFields     yes;
        executeControl  writeTime;
        writeControl    writeTime;
        patches         (wall);
    }

    yPlus
    {
        type            yPlus;
        libs            (fieldFunctionObjects);
        executeControl  writeTime;
        writeFields     yes;
        writeControl    writeTime;
    }
}
EOL

# Create fvSchemes
fv_schemes_file="$case_dir/system/fvSchemes"
cat > "$fv_schemes_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss upwind;
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
EOL

# Create fvSolution
fv_solution_file="$case_dir/system/fvSolution"
cat > "$fv_solution_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-6;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;

    residualControl
    {
        p               1e-4;
        U               1e-4;
        "(k|omega)"     1e-4;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        k               0.7;
        omega           0.7;
    }
}
EOL

# Create transportProperties
transport_file="$case_dir/constant/transportProperties"
cat > "$transport_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

transportModel  Newtonian;

nu              nu [0 2 -1 0 0 0 0] 3.3e-6;
EOL

# Create turbulenceProperties
mkdir -p "$case_dir/constant/turbulenceProperties"

# Create RASProperties (OpenFOAM 2.4.0 style)
turbulence_file="$case_dir/constant/turbulenceProperties/RASProperties"
cat > "$turbulence_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      RASProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

RASModel        kOmegaSST;

turbulence      on;

printCoeffs     on;
EOL

# Create initial conditions (p field)
p_file="$case_dir/0/p"
cat > "$p_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }

    wall
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# Create U (velocity) field
u_file="$case_dir/0/U"
cat > "$u_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0.4 0 0);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (0.4 0 0);
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            noSlip;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# Create k (turbulent kinetic energy) field
k_file="$case_dir/0/k"
cat > "$k_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.0006;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0.0006;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            kqRWallFunction;
        value           uniform 0.0006;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# Create omega (specific dissipation rate) field
omega_file="$case_dir/0/omega"
cat > "$omega_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform 100;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 100;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            omegaWallFunction;
        value           uniform 100;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# Create nut (turbulent viscosity) field
nut_file="$case_dir/0/nut"
cat > "$nut_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
    }

    wall
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# Create a simplified blockMeshDict with hardcoded values
echo "Creating simplified blockMeshDict with hard-coded values..."
block_mesh_file="$case_dir/system/blockMeshDict"
cat > "$block_mesh_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
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

// Convert mm to m
convertToMeters 0.001;

// Hard-coded values for 50% stenosis
vertices
(
    // Bottom face, inlet to stenosis start
    (0        -12.5     0)          // 0
    (65       -12.5     0)          // 1
    (65       -12.5     1)          // 2
    (0        -12.5     1)          // 3
    
    // Bottom face at stenosis
    (75       -6.25     0)          // 4
    (75       -6.25     1)          // 5
    
    // Bottom face, stenosis to outlet
    (85       -12.5     0)          // 6
    (85       -12.5     1)          // 7
    (150      -12.5     0)          // 8
    (150      -12.5     1)          // 9
    
    // Top face, inlet to stenosis start
    (0        12.5      0)          // 10
    (65       12.5      0)          // 11
    (65       12.5      1)          // 12
    (0        12.5      1)          // 13
    
    // Top face at stenosis
    (75       6.25      0)          // 14
    (75       6.25      1)          // 15
    
    // Top face, stenosis to outlet
    (85       12.5      0)          // 16
    (85       12.5      1)          // 17
    (150      12.5      0)          // 18
    (150      12.5      1)          // 19
);

blocks
(
    // Block 0: Inlet to stenosis start
    hex (0 1 11 10 3 2 12 13)
    (30 10 1)
    simpleGrading (1 1 1)
    
    // Block 1: Stenosis section (first half)
    hex (1 4 14 11 2 5 15 12)
    (20 10 1)
    simpleGrading (1 1 1)
    
    // Block 2: Stenosis section (second half)
    hex (4 6 16 14 5 7 17 15)
    (20 10 1)
    simpleGrading (1 1 1)
    
    // Block 3: Outlet section
    hex (6 8 18 16 7 9 19 17)
    (30 10 1)
    simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 10 13 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (8 9 19 18)
        );
    }
    wall
    {
        type wall;
        faces
        (
            // Bottom wall
            (0 1 2 3)
            (1 4 5 2)
            (4 6 7 5)
            (6 8 9 7)
            
            // Top wall
            (10 13 12 11)
            (11 12 15 14)
            (14 15 17 16)
            (16 17 19 18)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            // Front face (z=0)
            (0 1 11 10)
            (1 4 14 11)
            (4 6 16 14)
            (6 8 18 16)
            
            // Back face (z=thickness)
            (3 13 12 2)
            (2 12 15 5)
            (5 15 17 7)
            (7 17 19 9)
        );
    }
);

mergePatchPairs
(
);
EOL

# Create a case.foam file for ParaView
touch "$case_dir/case.foam"
echo "Created case.foam file for ParaView"

# Generate mesh
echo "Generating mesh with blockMesh..."
cd "$case_dir"
blockMesh > blockMesh.log 2>&1

if [ $? -ne 0 ]; then
    echo "blockMesh failed. Check log for details: blockMesh.log"
    cat blockMesh.log
    cd ..
    exit 1
fi
echo "blockMesh completed successfully"

# Check mesh quality
echo "Checking mesh quality..."
checkMesh > checkMesh.log 2>&1

# Run the simulation
echo "Starting OpenFOAM simulation..."
simpleFoam > simpleFoam.log 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo "Simulation completed successfully!"
    
    # Get the latest time directory
    latest_time=$(find . -maxdepth 1 -name "[0-9]*" | sort -n | tail -1)
    latest_time_value=$(basename "$latest_time")
    
    if [ ! -z "$latest_time_value" ]; then
        echo "Latest time directory: $latest_time_value"
        
        # Try to run some basic post-processing
        postProcess -func wallShearStress -time "$latest_time_value" > postWallShearStress.log 2>&1
        
        cd ..
        echo "Results available in directory: $case_dir"
    else
        cd ..
        echo "No time directories found. Simulation may have failed."
    fi
else
    echo "Simulation encountered errors. Check log file: simpleFoam.log"
    echo "Last 20 lines of the log file:"
    tail -n 20 simpleFoam.log
    cd ..
fi

echo "Job completed at $(date)"