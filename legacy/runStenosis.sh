#!/bin/bash
# runStenosis.sh - Simulate blood flow through an aorta with configurable stenosis level

# Function to display colorful messages
function echo_success() { echo -e "\033[0;32m✓ $1\033[0m"; }
function echo_info() { echo -e "\033[0;34mℹ $1\033[0m"; }
function echo_warning() { echo -e "\033[0;33m⚠ $1\033[0m"; }
function echo_error() { echo -e "\033[0;31m❌ $1\033[0m"; }

# Function to verify file creation
function verify_file() {
    if [ -f "$1" ]; then
        echo_success "Created file: $1"
    else
        echo_error "Failed to create file: $1"
        exit 1
    fi
}

# Source environment setup script
if [ -f "./setup_environment.sh" ]; then
    source ./setup_environment.sh
else
    echo_error "setup_environment.sh not found. Please create this file first."
    exit 1
fi

# Default values
stenosis_level=0.5
case_name="moderate"
mesh_size=0.5
max_iterations=2000
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
            echo "  --mesh SIZE       Mesh size for Gmsh, default: 0.5"
            echo "  --iterations NUM  Number of iterations, default: 2000"
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

# Calculate cell counts
TOTAL_CELLS_LENGTH=100
CELLS_FIRST=$(echo "($FIRST_LENGTH * $TOTAL_CELLS_LENGTH)/1" | bc)
CELLS_STENOSIS=$(echo "($STENOSIS_SECTION * $TOTAL_CELLS_LENGTH)/1" | bc)
CELLS_LAST=$(echo "($LAST_LENGTH * $TOTAL_CELLS_LENGTH)/1" | bc)
CELLS_CHECK=$(echo "$CELLS_FIRST + $CELLS_STENOSIS + $CELLS_LAST" | bc)

# Adjust if rounding created too few/many cells
if [ "$CELLS_CHECK" -ne "$TOTAL_CELLS_LENGTH" ]; then
    CELLS_FIRST=$(echo "$CELLS_FIRST + ($TOTAL_CELLS_LENGTH - $CELLS_CHECK)" | bc)
fi

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
    hex (0 1 11 10 3 2 12 13)
    ($CELLS_FIRST 15 1)
    simpleGrading (1 1 1)
    
    // Block 1: Stenosis section (first half)
    hex (1 4 14 11 2 5 15 12)
    ($CELLS_STENOSIS 15 1)
    simpleGrading (1 1 1)
    
    // Block 2: Stenosis section (second half)
    hex (4 6 16 14 5 7 17 15)
    ($CELLS_STENOSIS 15 1)
    simpleGrading (1 1 1)
    
    // Block 3: Outlet section
    hex (6 8 18 16 7 9 19 17)
    ($CELLS_LAST 15 1)
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

// ************************************************************************* //
EOL

verify_file "$block_mesh_file"

# Create decomposeParDict for parallel execution (optional)
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

# Create a Gmsh geometry file for the aorta with stenosis (still include this for reference)
echo_info "Creating Gmsh geometry file..."
gmsh_geo_file="$case_dir/aorta.geo"
cat > "$gmsh_geo_file" << EOL
// Parameters for aorta with atherosclerosis
diameter = 25.0;  // Aorta diameter in mm (average adult aorta ~25mm)
length = 150.0;   // Length of segment in mm
stenosis_level = $stenosis_level;  // 0.0 = no stenosis, 1.0 = complete occlusion
stenosis_length = 20.0;  // Length of the stenotic region in mm
stenosis_position = length/2;  // Position of stenosis center from inlet
mesh_size = $mesh_size;  // Default mesh size (smaller values create finer mesh)

// Calculate stenosis height based on stenosis level
stenosis_height = diameter * stenosis_level / 2;

// Points for upper wall
Point(1) = {0, diameter/2, 0, mesh_size};  // Inlet top
Point(2) = {stenosis_position - stenosis_length/2, diameter/2, 0, mesh_size};  // Start of stenosis top
Point(3) = {stenosis_position, diameter/2 - stenosis_height, 0, mesh_size/4};  // Peak of stenosis top (finer mesh)
Point(4) = {stenosis_position + stenosis_length/2, diameter/2, 0, mesh_size};  // End of stenosis top
Point(5) = {length, diameter/2, 0, mesh_size};  // Outlet top

// Points for lower wall
Point(6) = {0, -diameter/2, 0, mesh_size};  // Inlet bottom
Point(7) = {stenosis_position - stenosis_length/2, -diameter/2, 0, mesh_size};  // Start of stenosis bottom
Point(8) = {stenosis_position, -diameter/2 + stenosis_height, 0, mesh_size/4};  // Peak of stenosis bottom (finer mesh)
Point(9) = {stenosis_position + stenosis_length/2, -diameter/2, 0, mesh_size};  // End of stenosis bottom
Point(10) = {length, -diameter/2, 0, mesh_size};  // Outlet bottom

// Create splines for the walls (smoother than straight lines)
Spline(1) = {1, 2, 3, 4, 5};  // Upper wall
Spline(2) = {6, 7, 8, 9, 10};  // Lower wall

// Create inlet and outlet lines
Line(3) = {1, 6};  // Inlet
Line(4) = {5, 10}; // Outlet

// Create surface for the fluid domain
Line Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};

// Define physical groups for OpenFOAM boundary conditions
Physical Curve("inlet") = {3};
Physical Curve("outlet") = {4};
Physical Curve("wall") = {1, 2};
Physical Surface("fluid") = {1};

// Mesh settings for better quality
Mesh.RecombineAll = 1;  // Generate quadrilateral elements where possible
Mesh.Smoothing = 20;    // Mesh smoothing steps
Mesh.Algorithm = 6;     // Frontal-Delaunay for quads
EOL

verify_file "$gmsh_geo_file"

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

# Add function object for wall shear stress to controlDict
echo_info "Updating controlDict for wall shear stress calculation..."
wall_shear_dict="$case_dir/system/wallShearStressDict"
cat > "$wall_shear_dict" << EOL
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
    object      wallShearStressDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type            wallShearStress;
libs            (fieldFunctionObjects);
patches         (wall);
writeControl    writeTime;
writeFields     yes;
executeControl  writeTime;

// ************************************************************************* //
EOL

verify_file "$wall_shear_dict"

# STEP 4: MESH GENERATION AND SIMULATION SETUP
echo_info "Setting up mesh and simulation..."

# Option 1: Use blockMesh (more reliable)
echo_info "Generating mesh with blockMesh..."
blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1

if [ $? -ne 0 ]; then
    echo_error "blockMesh failed. Check log for details: $case_dir/blockMesh.log"
    cat "$case_dir/blockMesh.log"
    exit 1
fi
echo_success "blockMesh completed successfully"

# Option 2: Alternative to use Gmsh if blockMesh fails
if [ ! -d "$case_dir/constant/polyMesh" ] || [ ! -f "$case_dir/constant/polyMesh/points" ]; then
    echo_warning "blockMesh output not found, trying with Gmsh instead..."
    
    # Run Gmsh to generate the mesh
    echo_info "Generating mesh with Gmsh..."
    gmsh -2 "$case_dir/aorta.geo" -o "$case_dir/aorta.msh" > "$case_dir/gmsh.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo_error "Gmsh failed to generate the mesh. Check log: $case_dir/gmsh.log"
        cat "$case_dir/gmsh.log"
        exit 1
    fi
    
    echo_info "Converting Gmsh mesh to OpenFOAM format..."
    gmshToFoam "$case_dir/aorta.msh" -case "$case_dir" > "$case_dir/gmshToFoam.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo_error "gmshToFoam conversion failed. Check log: $case_dir/gmshToFoam.log"
        cat "$case_dir/gmshToFoam.log"
        exit 1
    fi
    
    # Update boundary file if needed
    echo_info "Updating boundary conditions after Gmsh conversion..."
    createPatch -overwrite -case "$case_dir" > "$case_dir/createPatch.log" 2>&1
fi

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

# STEP 5: RUN THE SIMULATION
echo_info "Starting OpenFOAM simulation for stenosis level $stenosis_level..."
echo_info "This may take several minutes depending on your computer..."

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
        
        # Create basic plots directory
        mkdir -p "plots"
        echo_info "Creating basic result plots in 'plots' directory..."
        
        if [ -f "plot_results.py" ]; then
            python3 plot_results.py --dirs "$case_dir" --all
            echo_success "Created visualization plots in 'plots' directory"
        else
            echo_warning "plot_results.py not found, skipping automatic plot generation"
        fi
    else
        echo_warning "No time directories found. Post-processing skipped."
    fi
    
    echo_info "You can visualize results with: paraFoam -case $case_dir"
else
    echo_error "Simulation encountered errors. Check log file: $case_dir/simpleFoam.log"
    echo_info "Last 20 lines of the log file:"
    tail -n 20 "$case_dir/simpleFoam.log"
    
    # Try to identify specific errors
    if grep -q "attempt to read beyond EOF" "$case_dir/simpleFoam.log"; then
        echo_error "EOF error detected in boundary conditions. This usually indicates a problem with the field files in the 0/ directory."
    fi
    if grep -q "could not find file" "$case_dir/simpleFoam.log"; then
        echo_error "Missing file error detected. Check that all required files exist."
    fi
    if grep -q "divergence" "$case_dir/simpleFoam.log"; then
        echo_error "Solution divergence detected. Try reducing relaxation factors in system/fvSolution."
    fi
fi

echo "======================================================="
echo "Simulation complete for stenosis level: $stenosis_level"
echo "======================================================="