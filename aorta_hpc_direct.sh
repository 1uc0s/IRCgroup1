#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -N aorta_stenosis
#PBS -j oe

# =========================================================
# DIRECT EXECUTION SCRIPT FOR AORTA CFD ON HPC
# This script skips runStenosisv3.sh and directly implements
# the necessary steps to run the simulation on HPC
# =========================================================

# Go to the directory from which the job was submitted
cd $PBS_O_WORKDIR
echo "Working directory: $(pwd)"

# Load required modules
echo "Loading required modules..."
module load tools/prod
module load mpi/intel-2018
module load openfoam/2.4.0

# Source OpenFOAM environment
echo "Setting up OpenFOAM environment..."
if [ -f "/apps/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc" ]; then
    echo "Using OpenFOAM at /apps/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc"
    source "/apps/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc"
elif [ -n "$FOAM_INST_DIR" ] && [ -f "$FOAM_INST_DIR/etc/bashrc" ]; then
    echo "Using FOAM_INST_DIR: $FOAM_INST_DIR/etc/bashrc"
    source "$FOAM_INST_DIR/etc/bashrc"
else
    echo "Could not find OpenFOAM environment file, but proceeding anyway..."
fi

# Verify OpenFOAM is available
if command -v simpleFoam &> /dev/null; then
    echo "OpenFOAM successfully loaded: $(which simpleFoam)"
    echo "OpenFOAM version: $(simpleFoam --version 2>&1 | head -1)"
else
    echo "ERROR: simpleFoam command not found. OpenFOAM not properly configured."
    exit 1
fi

# Set default values for simulation parameters
STENOSIS_LEVEL=${LEVEL:-0.5}
CASE_NAME=${NAME:-"moderate"}
MESH_SIZE=${MESH:-0.5}
MAX_ITERATIONS=${ITER:-1000}

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
# NOTE: This is copied from runStenosisv3.sh but streamlined for HPC use

echo "Creating OpenFOAM configuration files..."

# Create the controlDict
control_dict_file="$case_dir/system/controlDict"
cat > "$control_dict_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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

endTime         $MAX_ITERATIONS;

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

    forces
    {
        type            forces;
        libs            (forces);
        executeControl  writeTime;
        writeControl    writeTime;
        patches         (wall);
        rho             rhoInf;
        rhoInf          1060;
        CofR            (0 0 0);
        log             yes;
    }
    
    sampleLines
    {
        type            sets;
        libs            (sampling);
        executeControl  writeTime;
        writeControl    writeTime;
        
        sets
        {
            centerline
            {
                type    uniform;
                axis    distance;
                start   (0 0 0.5);
                end     (150 0 0.5);
                nPoints 500;
            }
            
            lowerWall
            {
                type    uniform;
                axis    distance;
                start   (0 -12.5 0.5);
                end     (150 -12.5 0.5);
                nPoints 500;
            }
            
            upperWall
            {
                type    uniform;
                axis    distance;
                start   (0 12.5 0.5);
                end     (150 12.5 0.5);
                nPoints 500;
            }
        }
        
        fields          (p U wallShearStress);
        interpolationScheme cellPoint;
        setFormat       raw;
    }
}
EOL

# Create fvSchemes
fv_schemes_file="$case_dir/system/fvSchemes"
cat > "$fv_schemes_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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
    grad(U)         cellLimited Gauss linear 1;
}

divSchemes
{
    default          none;
    div(phi,U)       bounded Gauss linearUpwind grad(U);
    turbulence       bounded Gauss upwind;
    div(phi,k)       \$turbulence;
    div(phi,omega)   \$turbulence;
    div(phi,nut)     \$turbulence;
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

wallDist
{
    method          meshWave;
}
EOL

# Create fvSolution
fv_solution_file="$case_dir/system/fvSolution"
cat > "$fv_solution_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }

    "(U|k|omega|nut)"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
        nSweeps         1;
    }

    "(U|k|omega|nut)Final"
    {
        \$U;
        relTol          0;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;
    
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
        U               0.6;
        k               0.5;
        omega           0.5;
        ".*"            0.5;
    }
}
EOL

# Create transportProperties
transport_file="$case_dir/constant/transportProperties"
cat > "$transport_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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

// Blood properties
nu              [0 2 -1 0 0 0 0] 3.3e-6;  // Kinematic viscosity of blood
rho             [1 -3 0 0 0 0 0] 1060;     // Density of blood in kg/m³
EOL

# Create turbulenceProperties
turbulence_file="$case_dir/constant/turbulenceProperties"
cat > "$turbulence_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType      RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
EOL

# Create initial conditions (p field)
p_file="$case_dir/0/p"
cat > "$p_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
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

# Create geometry and blockMeshDict
echo "Creating blockMeshDict..."

# Calculate stenosis radius and other parameters
AORTA_RADIUS=12.5 # mm
STENOSIS_RADIUS=$(echo "$AORTA_RADIUS * (1.0 - $STENOSIS_LEVEL)" | bc -l)
AORTA_LENGTH=150.0 # mm
STENOSIS_LENGTH=20.0 # mm
STENOSIS_POS=$(echo "$AORTA_LENGTH / 2" | bc -l)
INLET_POS=$(echo "$STENOSIS_POS - $STENOSIS_LENGTH / 2" | bc -l)
OUTLET_POS=$(echo "$STENOSIS_POS + $STENOSIS_LENGTH / 2" | bc -l)

echo "Geometry parameters:"
echo "- Aorta radius: $AORTA_RADIUS mm"
echo "- Stenosis radius: $STENOSIS_RADIUS mm (${STENOSIS_LEVEL} occlusion)"
echo "- Stenosis position: $STENOSIS_POS mm"

# Calculate mesh refinement
MESH_REFINEMENT_FACTOR=1.0
MESH_REFINEMENT_FACTOR=$(echo "$MESH_REFINEMENT_FACTOR / $MESH_SIZE" | bc -l)

# Mesh cells calculation
BASE_CELLS_LENGTH=80  # Reduced for HPC
BASE_CELLS_CROSS=12  

# Calculate cells
TOTAL_CELLS_LENGTH=$(echo "($BASE_CELLS_LENGTH * $MESH_REFINEMENT_FACTOR)/1" | bc)
CELLS_CROSS=$(echo "($BASE_CELLS_CROSS * $MESH_REFINEMENT_FACTOR)/1" | bc)
if [ "$CELLS_CROSS" -lt 6 ]; then
    CELLS_CROSS=6
fi

# Calculate length proportions
FIRST_LENGTH=$(echo "$INLET_POS / $AORTA_LENGTH" | bc -l)
STENOSIS_SECTION=$(echo "$STENOSIS_LENGTH / $AORTA_LENGTH" | bc -l)
LAST_LENGTH=$(echo "1.0 - $FIRST_LENGTH - $STENOSIS_SECTION" | bc -l)

# Calculate cell counts for each section
CELLS_INLET=$(echo "($FIRST_LENGTH * $TOTAL_CELLS_LENGTH * 0.8)/1" | bc)
CELLS_STENOSIS=$(echo "($STENOSIS_SECTION * $TOTAL_CELLS_LENGTH * 1.5)/1" | bc)
CELLS_OUTLET=$(echo "($OUTLET_PROP * $TOTAL_CELLS_LENGTH * 0.8)/1" | bc)

# Verify cell distribution
CELLS_SUM=$(echo "$CELLS_INLET + $CELLS_STENOSIS + $CELLS_OUTLET" | bc)
if [ "$CELLS_SUM" -ne "$TOTAL_CELLS_LENGTH" ]; then
    CELLS_INLET=$(echo "$CELLS_INLET + ($TOTAL_CELLS_LENGTH - $CELLS_SUM)" | bc)
fi

# Create blockMeshDict file - simplified for HPC
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
scale 0.001;

vertices
(
    // Bottom face, inlet to stenosis start
    (0              -$AORTA_RADIUS      0)          // 0
    ($INLET_POS     -$AORTA_RADIUS      0)          // 1
    ($INLET_POS     -$AORTA_RADIUS      1)          // 2
    (0              -$AORTA_RADIUS      1)          // 3
    
    // Bottom face at stenosis
    ($STENOSIS_POS  -$STENOSIS_RADIUS   0)          // 4
    ($STENOSIS_POS  -$STENOSIS_RADIUS   1)          // 5
    
    // Bottom face, stenosis to outlet
    ($OUTLET_POS    -$AORTA_RADIUS      0)          // 6
    ($OUTLET_POS    -$AORTA_RADIUS      1)          // 7
    ($AORTA_LENGTH  -$AORTA_RADIUS      0)          // 8
    ($AORTA_LENGTH  -$AORTA_RADIUS      1)          // 9
    
    // Top face, inlet to stenosis start
    (0              $AORTA_RADIUS       0)          // 10
    ($INLET_POS     $AORTA_RADIUS       0)          // 11
    ($INLET_POS     $AORTA_RADIUS       1)          // 12
    (0              $AORTA_RADIUS       1)          // 13
    
    // Top face at stenosis
    ($STENOSIS_POS  $STENOSIS_RADIUS    0)          // 14
    ($STENOSIS_POS  $STENOSIS_RADIUS    1)          // 15
    
    // Top face, stenosis to outlet
    ($OUTLET_POS    $AORTA_RADIUS       0)          // 16
    ($OUTLET_POS    $AORTA_RADIUS       1)          // 17
    ($AORTA_LENGTH  $AORTA_RADIUS       0)          // 18
    ($AORTA_LENGTH  $AORTA_RADIUS       1)          // 19
);

blocks
(
    // Block 0: Inlet to stenosis start
    hex (0 1 11 10 3 2 12 13)
    ($CELLS_INLET $CELLS_CROSS 1)
    simpleGrading (1 1 1)
    
    // Block 1: Stenosis section (first half)
    hex (1 4 14 11 2 5 15 12)
    ($CELLS_STENOSIS $CELLS_CROSS 1)
    simpleGrading (1 1 1)
    
    // Block 2: Stenosis section (second half)
    hex (4 6 16 14 5 7 17 15)
    ($CELLS_STENOSIS $CELLS_CROSS 1)
    simpleGrading (1 1 1)
    
    // Block 3: Outlet section
    hex (6 8 18 16 7 9 19 17)
    ($CELLS_OUTLET $CELLS_CROSS 1)
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

# Create decomposeParDict for parallel execution
decompose_file="$case_dir/system/decomposeParDict"
cat > "$decompose_file" << EOL
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
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;  // Use 4 cores

method          scotch;

// ************************************************************************* //
EOL

# Create a case.foam file for ParaView
touch "$case_dir/case.foam"
echo "Created case.foam file for ParaView"

# Generate mesh
echo "Generating mesh with blockMesh..."
cd "$case_dir"
blockMesh > blockMesh.log 2>&1

if [ $? -ne 0 ]; then
    echo "blockMesh failed. Check log for details: $case_dir/blockMesh.log"
    cat blockMesh.log
    cd ..
    exit 1
fi
echo "blockMesh completed successfully"

# Check mesh quality
echo "Checking mesh quality..."
checkMesh > checkMesh.log 2>&1

# Start the simulation
echo "Starting OpenFOAM simulation for stenosis level $STENOSIS_LEVEL..."
echo "This may take several minutes depending on your computer..."

# Run simpleFoam
simpleFoam > simpleFoam.log 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo "Simulation completed successfully!"
    
    # Get the latest time directory
    latest_time=$(find . -maxdepth 1 -name "[0-9]*" | sort -n | tail -1)
    latest_time_value=$(basename "$latest_time")
    
    if [ ! -z "$latest_time_value" ]; then
        echo "Latest time directory: $latest_time_value"
        
        # Calculate wall shear stress
        postProcess -func wallShearStress -time "$latest_time_value" > postWallShearStress.log 2>&1
        
        # Sample data along centerline and walls
        postProcess -func sample -time "$latest_time_value" > postSample.log 2>&1
        
        cd ..
        echo "Results available in directory: $case_dir"
        
        # Create plots if Python is available
        # if command -v python3 &> /dev/null && [ -f "plot_results.py" ]; then
        #     mkdir -p "plots"
        #     python3 plot_results.py --dirs "$case_dir" --all
        # fi
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