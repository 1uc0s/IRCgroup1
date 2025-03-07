#!/bin/bash
# runStenosis_hpc.sh - Standalone script for running aorta stenosis simulation on HPC
# This version is stripped down to work on HPC environments and take command line parameters

# Function to display colorful messages (if terminal supports it)
function echo_success() { echo "✓ $1"; }
function echo_info() { echo "ℹ $1"; }
function echo_warning() { echo "⚠ $1"; }
function echo_error() { echo "❌ $1"; }

# Function to verify file creation
function verify_file() {
    if [ -f "$1" ]; then
        echo_success "Created file: $1"
    else
        echo_error "Failed to create file: $1"
        exit 1
    fi
}

# Default values
stenosis_level=0.5
case_name="moderate"
mesh_size=0.5
max_iterations=1000
write_interval=200

# Process command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --level)
            stenosis_level="$2"
            shift 2
            ;;
        --name)
            case_name="$2"
            shift 2
            ;;
        --mesh)
            mesh_size="$2"
            shift 2
            ;;
        --iterations)
            max_iterations="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --level VALUE     Stenosis level (0.0-0.9), default: 0.5"
            echo "  --name NAME       Name for the case, default: moderate"
            echo "  --mesh SIZE       Mesh size factor (smaller = finer), default: 0.5"
            echo "  --iterations NUM  Number of iterations, default: 1000"
            echo "  --help            Display this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$0 --help' for usage information."
            exit 1
            ;;
    esac
done

# Create case directory with explicit paths
case_dir="aorta_simulation_${case_name}"
echo_info "Creating simulation directory structure: $case_dir"

# Create directories explicitly to avoid brace expansion issues
mkdir -p "$case_dir"
mkdir -p "$case_dir/0"
mkdir -p "$case_dir/constant"
mkdir -p "$case_dir/system"

# Verify directories were created
for dir in "$case_dir" "$case_dir/0" "$case_dir/constant" "$case_dir/system"; do
    if [ -d "$dir" ]; then
        echo_success "Created directory: $dir"
    else
        echo_error "Failed to create directory: $dir"
        exit 1
    fi
done

echo "======================================================="
echo "Processing stenosis level: $case_name ($stenosis_level)"
echo "Mesh size factor: $mesh_size"
echo "======================================================="

# STEP 1: CREATE ALL OPENFOAM CONFIGURATION FILES FIRST

# Create system files
echo_info "Creating OpenFOAM configuration files..."

# Create the controlDict with wall shear stress function enabled
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

endTime         $max_iterations;  // Iterations for steady-state solver

deltaT          1;     // Time step for steady solver is iteration step

writeControl    timeStep;

writeInterval   $write_interval;   // Write results every $write_interval iterations

purgeWrite      3;     // Keep only the 3 most recent time directories

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   8;

runTimeModifiable true;

// Function objects for post-processing
functions
{
    // Calculate wall shear stress
    wallShearStress
    {
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        writeFields     yes;
        executeControl  writeTime;
        writeControl    writeTime;
        patches         (wall);
    }

    // Calculate y+ to check mesh quality near walls
    yPlus
    {
        type            yPlus;
        libs            (fieldFunctionObjects);
        executeControl  writeTime;
        writeFields     yes;
        writeControl    writeTime;
    }

    // Calculate forces on the wall
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

    // Write solver residuals
    residuals
    {
        type            solverInfo;
        libs            (utilityFunctionObjects);
        executeControl  timeStep;
        writeControl    timeStep;
        fields          (".*");
    }
    
    // Sample fields along lines
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

verify_file "$control_dict_file"

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

verify_file "$fv_schemes_file"

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
    
    // Adding reference pressure settings
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
        ".*"            0.5;
    }
}

// ************************************************************************* //
EOL

verify_file "$fv_solution_file"

# Create constant directory files
echo_info "Creating physical properties files..."

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

verify_file "$turbulence_file"

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

verify_file "$transport_file"

# STEP 2: CREATE INITIAL CONDITIONS (0 directory)
echo_info "Creating initial conditions..."

# p (pressure)
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

// ************************************************************************* //
EOL

verify_file "$p_file"

# U (velocity)
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

// ************************************************************************* //
EOL

verify_file "$u_file"

# k (turbulent kinetic energy)
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

// ************************************************************************* //
EOL

verify_file "$k_file"

# omega (specific dissipation rate)
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

// ************************************************************************* //
EOL

verify_file "$omega_file"

# nut (turbulent viscosity)
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

// ************************************************************************* //
EOL

verify_file "$nut_file"

# STEP 3: CREATE GEOMETRY AND MESH

# Create a fixed blockMeshDict file without using #calc or #codeStream
echo_info "Creating blockMesh dictionary..."
block_mesh_file="$case_dir/system/blockMeshDict"

# Calculate stenosis radius (safer to do it here in the script)
AORTA_RADIUS=12.5 # mm
STENOSIS_RADIUS=$(echo "$AORTA_RADIUS * (1.0 - $stenosis_level)" | bc -l)
AORTA_LENGTH=150.0 # mm
STENOSIS_LENGTH=20.0 # mm
STENOSIS_POS=$(echo "$AORTA_LENGTH / 2" | bc -l)
INLET_POS=$(echo "$STENOSIS_POS - $STENOSIS_LENGTH / 2" | bc -l)
OUTLET_POS=$(echo "$STENOSIS_POS + $STENOSIS_LENGTH / 2" | bc -l)

# Calculate proportions for cell distribution
FIRST_LENGTH=$(echo "$INLET_POS / $AORTA_LENGTH" | bc -l)
STENOSIS_SECTION=$(echo "$STENOSIS_LENGTH / $AORTA_LENGTH" | bc -l)
LAST_LENGTH=$(echo "1.0 - $FIRST_LENGTH - $STENOSIS_SECTION" | bc -l)

# Calculate mesh refinement parameters
MESH_REFINEMENT_FACTOR=$(echo "scale=2; 1.0/$mesh_size" | bc)
echo_info "Mesh refinement factor: $MESH_REFINEMENT_FACTOR"

# Calculate cells based on mesh refinement
BASE_CELLS_LENGTH=100
BASE_CELLS_CROSS=16

TOTAL_CELLS_LENGTH=$(echo "($BASE_CELLS_LENGTH * $MESH_REFINEMENT_FACTOR)/1" | bc)
CELLS_CROSS=$(echo "($BASE_CELLS_CROSS * $MESH_REFINEMENT_FACTOR)/1" | bc)
# Ensure at least 6 cells across
if [ "$CELLS_CROSS" -lt 6 ]; then
    CELLS_CROSS=6
fi

# Calculate cell counts for each section with refinement at stenosis
CELLS_INLET=$(echo "($FIRST_LENGTH * $TOTAL_CELLS_LENGTH * 0.8)/1" | bc)
CELLS_STENOSIS=$(echo "($STENOSIS_PROP * $TOTAL_CELLS_LENGTH * 1.5)/1" | bc)
CELLS_OUTLET=$(echo "($OUTLET_PROP * $TOTAL_CELLS_LENGTH * 0.8)/1" | bc)

# Verify cell distribution adds up to total
CELLS_SUM=$(echo "$CELLS_INLET + $CELLS_STENOSIS + $CELLS_OUTLET" | bc)
if [ "$CELLS_SUM" -ne "$TOTAL_CELLS_LENGTH" ]; then
    CELLS_INLET=$(echo "$CELLS_INLET + ($TOTAL_CELLS_LENGTH - $CELLS_SUM)" | bc)
fi

echo_info "Cell counts - Length: $TOTAL_CELLS_LENGTH (inlet: $CELLS_INLET, stenosis: $CELLS_STENOSIS, outlet: $CELLS_OUTLET), Cross: $CELLS_CROSS"

# Create the blockMeshDict with explicitly calculated values
cat > "$block_mesh_file" << EOL
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
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convert mm to m
scale 0.001;

// Stenosis parameters (pre-calculated)
// Radius: $AORTA_RADIUS mm
// Stenosis level: $stenosis_level
// Stenosis radius: $STENOSIS_RADIUS mm

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
    // Grading: 1/2 in x-direction to refine cells toward stenosis
    // Grading: 4 in y-direction for wall refinement
    hex (0 1 11 10 3 2 12 13)
    ($CELLS_INLET $CELLS_CROSS 1)
    simpleGrading (2 4 1)
    
    // Block 1: Stenosis section (first half)
    // Symmetric grading to refine at stenosis peak
    hex (1 4 14 11 2 5 15 12)
    ($CELLS_STENOSIS $CELLS_CROSS 1)
    simpleGrading (0.5 4 1)
    
    // Block 2: Stenosis section (second half)
    // Symmetric grading to refine at stenosis peak
    hex (4 6 16 14 5 7 17 15)
    ($CELLS_STENOSIS $CELLS_CROSS 1)
    simpleGrading (2 4 1)
    
    // Block 3: Outlet section
    // Grading: 2 in x-direction to refine cells near stenosis
    hex (6 8 18 16 7 9 19 17)
    ($CELLS_OUTLET $CELLS_CROSS 1)
    simpleGrading (0.5 4 1)
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

// ************************************************************************* //
EOL

verify_file "$block_mesh_file"

# Create decomposeParDict for parallel execution
decompose_file="$case_dir/system/decomposeParDict"
cat > "$decompose_file" << EOL
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
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;  // Adjust based on available cores

method          scotch;  // No geometric info needed

// ************************************************************************* //
EOL

verify_file "$decompose_file"

# Create sampleDict for post-processing
echo_info "Creating sampling dictionary for centerline data extraction..."
sample_dict_file="$case_dir/system/sampleDict"
cat > "$sample_dict_file" << EOL
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
    object      sampleDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type            sets;
libs            (sampling);
interpolationScheme cellPoint;
setFormat       raw;

fields          ( p U wallShearStress );

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

// ************************************************************************* //
EOL

verify_file "$sample_dict_file"

# STEP 4: MESH GENERATION AND SIMULATION SETUP
echo_info "Setting up mesh and simulation..."

# Generate mesh with blockMesh
echo_info "Generating mesh with blockMesh..."
blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1

if [ $? -ne 0 ]; then
    echo_error "blockMesh failed. Check log for details: $case_dir/blockMesh.log"
    cat "$case_dir/blockMesh.log"
    exit 1
fi
echo_success "blockMesh completed successfully"

# Check mesh quality
echo_info "Checking mesh quality..."
checkMesh -case "$case_dir" > "$case_dir/checkMesh.log" 2>&1

# Check for mesh errors
if grep -q "FAILED" "$case_dir/checkMesh.log"; then
    echo_warning "Mesh check found issues. See details in: $case_dir/checkMesh.log"
    echo_info "This may impact simulation accuracy, but we'll proceed anyway."
else
    echo_success "Mesh check passed successfully"
fi

# Create a case.foam file for ParaView
touch "$case_dir/case.foam"
echo_success "Created case.foam file for ParaView"

# STEP 5: RUN THE SIMULATION
echo_info "Starting OpenFOAM simulation for stenosis level $stenosis_level..."
echo_info "Running for $max_iterations iterations..."

# Run in serial mode
simpleFoam -case "$case_dir" > "$case_dir/simpleFoam.log" 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo_success "Simulation completed successfully!"
    echo_info "Results available in directory: $case_dir"
    
    # Run post-processing for wall shear stress
    echo_info "Running post-processing for visualization data..."
    
    # Get the latest time directory
    latest_time=$(find "$case_dir" -maxdepth 1 -name "[0-9]*" | sort -n | tail -1)
    latest_time_value=$(basename "$latest_time")
    
    if [ ! -z "$latest_time_value" ]; then
        echo_info "Using latest time directory: $latest_time_value"
        
        # Calculate wall shear stress using the wallShearStress function object
        postProcess -func wallShearStress -case "$case_dir" -time "$latest_time_value" > "$case_dir/postWallShearStress.log" 2>&1
        
        # Sample data along centerline and walls
        postProcess -func sample -case "$case_dir" -time "$latest_time_value" > "$case_dir/postSample.log" 2>&1
    else
        echo_warning "No time directories found. Post-processing skipped."
    fi
else
    echo_error "Simulation encountered errors. Check log file: $case_dir/simpleFoam.log"
    echo_info "Last 20 lines of the log file:"
    tail -n 20 "$case_dir/simpleFoam.log"
    exit 1
fi

echo "======================================================="
echo "Simulation complete for stenosis level: $stenosis_level"
echo "======================================================="

# Calculate and output metrics for the convergence study
# This is particularly helpful for the mesh convergence study

# Extract pressure drop
echo_info "Calculating key metrics for convergence study..."

# Get centerline data for pressure
CENTERLINE_FILE="${latest_time}/centerline_p_raw.xy"
if [ -f "$case_dir/$CENTERLINE_FILE" ]; then
    # Calculate pressure drop (inlet - outlet)
    INLET_PRESSURE=$(head -n 10 "$case_dir/$CENTERLINE_FILE" | awk '{sum+=$2; count++} END {print sum/count}')
    OUTLET_PRESSURE=$(tail -n 10 "$case_dir/$CENTERLINE_FILE" | awk '{sum+=$2; count++} END {print sum/count}')
    PRESSURE_DROP=$(echo "$INLET_PRESSURE - $OUTLET_PRESSURE" | bc -l)
    
    echo "KEY METRIC: Pressure Drop = $PRESSURE_DROP Pa"
    echo "$PRESSURE_DROP" > "$case_dir/pressure_drop.txt"
else
    echo_warning "Pressure data not found, skipping pressure drop calculation"
fi

# Get wall shear stress data
WSS_FILE="${latest_time}/lowerWall_wallShearStress_raw.xy"
if [ -f "$case_dir/$WSS_FILE" ]; then
    # Calculate maximum wall shear stress
    MAX_WSS=$(awk 'BEGIN{max=0} {mag=sqrt($2*$2+$3*$3+$4*$4); if(mag>max) max=mag} END{print max}' "$case_dir/$WSS_FILE")
    
    echo "KEY METRIC: Maximum Wall Shear Stress = $MAX_WSS Pa"
    echo "$MAX_WSS" > "$case_dir/max_wss.txt"
else 
    echo_warning "Wall shear stress data not found, skipping WSS calculation"
fi

# Output mesh size for reference
CELL_COUNT=$(grep "cells:" "$case_dir/checkMesh.log" | awk '{print $2}')
echo "Mesh Size: $mesh_size (Cell count: $CELL_COUNT)"
echo "$mesh_size $CELL_COUNT" > "$case_dir/mesh_info.txt"

echo "Convergence study metrics written to: $case_dir/pressure_drop.txt and $case_dir/max_wss.txt"