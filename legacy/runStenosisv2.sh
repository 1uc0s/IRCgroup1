#!/bin/bash
# runStenosis.sh - Unified script for aorta stenosis simulation with blockMesh and Gmsh support
# Features improved mesh control, field consistency, and robust error handling

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

# Source environment setup script if available
if [ -f "./setup_environment.sh" ]; then
    source ./setup_environment.sh
    echo_success "Sourced OpenFOAM environment"
else
    echo_warning "setup_environment.sh not found. Proceeding with current environment."
    # Check if OpenFOAM is available
    if ! command -v simpleFoam &> /dev/null; then
        echo_error "OpenFOAM not found in PATH. Make sure OpenFOAM is installed and sourced."
        exit 1
    fi
fi

# Default values
stenosis_level=0.5
case_name="moderate"
mesh_size=0.5
max_iterations=2000
write_interval=200
mesher="block"  # Options: "block" or "gmsh"
mesh_quality="standard" # Options: "standard", "fine", or "coarse"

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
        --mesher)
            mesher="$2"
            if [[ "$mesher" != "block" && "$mesher" != "gmsh" ]]; then
                echo_error "Invalid mesher option. Use 'block' or 'gmsh'."
                exit 1
            fi
            shift 2
            ;;
        --quality)
            mesh_quality="$2"
            if [[ "$mesh_quality" != "coarse" && "$mesh_quality" != "standard" && "$mesh_quality" != "fine" ]]; then
                echo_error "Invalid mesh quality. Use 'coarse', 'standard', or 'fine'."
                exit 1
            fi
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
            echo "  --mesh SIZE       Base mesh size, default: 0.5 (smaller = finer)"
            echo "  --mesher TYPE     Meshing method: 'block' (default) or 'gmsh'"
            echo "  --quality LEVEL   Mesh quality: 'coarse', 'standard' (default), or 'fine'"
            echo "  --iterations NUM  Number of solver iterations, default: 2000"
            echo "  --help            Display this help message"
            exit 0
            ;;
        *)
            echo_error "Unknown option: $1"
            echo "Run '$0 --help' for usage information."
            exit 1
            ;;
    esac
done

# Create case directory with explicit paths
case_dir="aorta_simulation_${case_name}"
echo_info "Creating simulation in: $case_dir"

# Create clean directory structure
if [ -d "$case_dir" ]; then
    echo_warning "Directory $case_dir already exists. Cleaning up..."
    rm -rf "$case_dir"
fi

# Create directories explicitly
mkdir -p "$case_dir"
mkdir -p "$case_dir/0"
mkdir -p "$case_dir/constant"
mkdir -p "$case_dir/system"
mkdir -p "$case_dir/constant/polyMesh"

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
echo "Using mesher: $mesher with $mesh_quality quality"
echo "======================================================="

# STEP 1: CREATE ALL OPENFOAM CONFIGURATION FILES 

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

# Create fvSolution - adjust relaxation factors for better stability
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
        U               0.6;  // Slightly reduced for stability
        k               0.5;
        omega           0.5;
        ".*"            0.5;
    }
}

// ************************************************************************* //
EOL

verify_file "$fv_solution_file"

# Create transportProperties with blood properties
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

# Create turbulenceProperties with kOmegaSST model
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

# STEP 3: CREATE GEOMETRIC PARAMETERS AND MESH

# Calculate geometry parameters
AORTA_RADIUS=12.5 # mm
STENOSIS_RADIUS=$(echo "$AORTA_RADIUS * (1.0 - $stenosis_level)" | bc -l)
AORTA_LENGTH=150.0 # mm
STENOSIS_LENGTH=20.0 # mm
STENOSIS_POS=$(echo "$AORTA_LENGTH / 2" | bc -l)
INLET_POS=$(echo "$STENOSIS_POS - $STENOSIS_LENGTH / 2" | bc -l)
OUTLET_POS=$(echo "$STENOSIS_POS + $STENOSIS_LENGTH / 2" | bc -l)

echo_info "Geometry parameters:"
echo_info "- Aorta radius: $AORTA_RADIUS mm"
echo_info "- Stenosis radius: $STENOSIS_RADIUS mm (${stenosis_level} occlusion)"
echo_info "- Stenosis position: $STENOSIS_POS mm"

# Set mesh refinement parameters based on quality
if [ "$mesh_quality" == "fine" ]; then
    MESH_REFINEMENT_FACTOR=2.0
    STENOSIS_GRADING=0.25  # Finer mesh at stenosis
elif [ "$mesh_quality" == "coarse" ]; then
    MESH_REFINEMENT_FACTOR=0.5
    STENOSIS_GRADING=0.5   # Less refinement
else  # standard
    MESH_REFINEMENT_FACTOR=1.0
    STENOSIS_GRADING=0.3   # Balanced
fi

# Scale with user-provided mesh size
MESH_REFINEMENT_FACTOR=$(echo "$MESH_REFINEMENT_FACTOR / $mesh_size" | bc -l)

# BlockMesh parameters
BASE_CELLS_LENGTH=100
BASE_CELLS_CROSS=16
CELLS_DEPTH=1

# Calculate cells based on mesh refinement
TOTAL_CELLS_LENGTH=$(echo "($BASE_CELLS_LENGTH * $MESH_REFINEMENT_FACTOR)/1" | bc)
CELLS_CROSS=$(echo "($BASE_CELLS_CROSS * $MESH_REFINEMENT_FACTOR)/1" | bc)
# Ensure at least 6 cells across radius
if [ "$CELLS_CROSS" -lt 6 ]; then
    CELLS_CROSS=6
fi

# Calculate length sections for better cell distribution around stenosis
FIRST_PROP=$(echo "$INLET_POS / $AORTA_LENGTH" | bc -l)
STENOSIS_PROP=$(echo "$STENOSIS_LENGTH / $AORTA_LENGTH" | bc -l)
OUTLET_PROP=$(echo "1.0 - $FIRST_PROP - $STENOSIS_PROP" | bc -l)

# Calculate cell counts for each section with refinement at stenosis
CELLS_INLET=$(echo "($FIRST_PROP * $TOTAL_CELLS_LENGTH * 0.8)/1" | bc)
CELLS_STENOSIS=$(echo "($STENOSIS_PROP * $TOTAL_CELLS_LENGTH * 1.5)/1" | bc)
CELLS_OUTLET=$(echo "($OUTLET_PROP * $TOTAL_CELLS_LENGTH * 0.8)/1" | bc)

# Verify cell distribution adds up to total
CELLS_SUM=$(echo "$CELLS_INLET + $CELLS_STENOSIS + $CELLS_OUTLET" | bc)
if [ "$CELLS_SUM" -ne "$TOTAL_CELLS_LENGTH" ]; then
    CELLS_INLET=$(echo "$CELLS_INLET + ($TOTAL_CELLS_LENGTH - $CELLS_SUM)" | bc)
fi

echo_info "Mesh parameters:"
echo_info "- Refinement factor: $MESH_REFINEMENT_FACTOR"
echo_info "- Total axial cells: $TOTAL_CELLS_LENGTH (inlet: $CELLS_INLET, stenosis: $CELLS_STENOSIS, outlet: $CELLS_OUTLET)"
echo_info "- Radial cells: $CELLS_CROSS"

# Create a blockMeshDict with grading to refine near stenosis and walls
if [ "$mesher" == "block" ]; then
    echo_info "Creating blockMesh definition with grading..."
    block_mesh_file="$case_dir/system/blockMeshDict"
    
    # Calculate grading parameters for stenosis refinement
    # This creates finer cells at the stenosis
    Y_GRADING=4  # Radial grading for finer cells near wall
    
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

// Aorta parameters
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
    // Grading: $Y_GRADING in y-direction for wall refinement
    hex (0 1 11 10 3 2 12 13)
    ($CELLS_INLET $CELLS_CROSS 1)
    simpleGrading (2 $Y_GRADING 1)
    
    // Block 1: Stenosis section (first half)
    // Symmetric grading to refine at stenosis peak
    hex (1 4 14 11 2 5 15 12)
    ($CELLS_STENOSIS $CELLS_CROSS 1)
    simpleGrading (0.5 $Y_GRADING 1)
    
    // Block 2: Stenosis section (second half)
    // Symmetric grading to refine at stenosis peak
    hex (4 6 16 14 5 7 17 15)
    ($CELLS_STENOSIS $CELLS_CROSS 1)
    simpleGrading (2 $Y_GRADING 1)
    
    // Block 3: Outlet section
    // Grading: 2 in x-direction to refine cells near stenosis
    hex (6 8 18 16 7 9 19 17)
    ($CELLS_OUTLET $CELLS_CROSS 1)
    simpleGrading (0.5 $Y_GRADING 1)
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
fi

# Create Gmsh geometry with variable mesh refinement
echo_info "Creating Gmsh geometry with variable refinement..."
gmsh_geo_file="$case_dir/aorta.geo"

# Adjust mesh parameters for Gmsh
if [ "$mesh_quality" == "fine" ]; then
    GMSH_BASE_SIZE=0.3
    GMSH_STENOSIS_SIZE=0.1
elif [ "$mesh_quality" == "coarse" ]; then
    GMSH_BASE_SIZE=1.0
    GMSH_STENOSIS_SIZE=0.4
else  # standard
    GMSH_BASE_SIZE=0.6
    GMSH_STENOSIS_SIZE=0.2
fi

# Scale with user-provided mesh size
GMSH_BASE_SIZE=$(echo "$GMSH_BASE_SIZE * $mesh_size" | bc -l)
GMSH_STENOSIS_SIZE=$(echo "$GMSH_STENOSIS_SIZE * $mesh_size" | bc -l)

cat > "$gmsh_geo_file" << EOL
// Parameters for aorta with atherosclerosis
diameter = 25.0;  // Aorta diameter in mm (average adult aorta ~25mm)
length = 150.0;   // Length of segment in mm
stenosis_level = $stenosis_level;  // 0.0 = no stenosis, 1.0 = complete occlusion
stenosis_length = 20.0;  // Length of the stenotic region in mm
stenosis_position = length/2;  // Position of stenosis center from inlet
base_mesh_size = $GMSH_BASE_SIZE;  // Base mesh size (smaller values create finer mesh)
stenosis_mesh_size = $GMSH_STENOSIS_SIZE; // Fine mesh at stenosis

// Calculate stenosis height based on stenosis level
stenosis_height = diameter * stenosis_level / 2;

// Points for upper wall
Point(1) = {0, diameter/2, 0, base_mesh_size};  // Inlet top
Point(2) = {stenosis_position - stenosis_length/2, diameter/2, 0, base_mesh_size};  // Start of stenosis top
Point(3) = {stenosis_position, diameter/2 - stenosis_height, 0, stenosis_mesh_size};  // Peak of stenosis top (finer mesh)
Point(4) = {stenosis_position + stenosis_length/2, diameter/2, 0, base_mesh_size};  // End of stenosis top
Point(5) = {length, diameter/2, 0, base_mesh_size};  // Outlet top

// Points for lower wall
Point(6) = {0, -diameter/2, 0, base_mesh_size};  // Inlet bottom
Point(7) = {stenosis_position - stenosis_length/2, -diameter/2, 0, base_mesh_size};  // Start of stenosis bottom
Point(8) = {stenosis_position, -diameter/2 + stenosis_height, 0, stenosis_mesh_size};  // Peak of stenosis bottom (finer mesh)
Point(9) = {stenosis_position + stenosis_length/2, -diameter/2, 0, base_mesh_size};  // End of stenosis bottom
Point(10) = {length, -diameter/2, 0, base_mesh_size};  // Outlet bottom

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

// Create mesh size field to refine near stenosis
Field[1] = Distance;
Field[1].NodesList = {3, 8}; // Points at stenosis peak

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = stenosis_mesh_size;
Field[2].LcMax = base_mesh_size;
Field[2].DistMin = 0;
Field[2].DistMax = stenosis_length;

// Use the minimum of all fields as the mesh size
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;

// Mesh settings for better quality
Mesh.Algorithm = 6;     // Frontal-Delaunay for quads (more robust)
Mesh.RecombineAll = 1;  // Generate quadrilateral elements where possible
Mesh.Smoothing = 20;    // Mesh smoothing steps
Mesh.ElementOrder = 1;  // First-order elements
Mesh.SecondOrderLinear = 0; // Linear interpolation for second-order elements
Mesh.CharacteristicLengthExtendFromBoundary = 1;
Mesh.CharacteristicLengthMin = stenosis_mesh_size / 2;
Mesh.CharacteristicLengthMax = base_mesh_size * 1.5;
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

# Create createPatchDict for boundary conversion if needed
createpatch_file="$case_dir/system/createPatchDict"
cat > "$createpatch_file" << EOL
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
    object      createPatchDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Do a synchronisation of coupled points after creation of any patches.
// Note: this does not work with points that are on multiple coupled patches
//       with transformations (i.e. cyclics).
pointSync false;

// Patches to create.
patches
(
    {
        // Name of new patch
        name inlet;

        // Dictionary to construct new patch from
        patchInfo
        {
            type patch;
        }

        // How to construct: either from 'patches' or 'set'
        constructFrom patches;

        // If constructFrom = patches : names of patches. Wildcards allowed.
        patches (inlet);
    }

    {
        // Name of new patch
        name outlet;

        // Dictionary to construct new patch from
        patchInfo
        {
            type patch;
        }

        // How to construct: either from 'patches' or 'set'
        constructFrom patches;

        // If constructFrom = patches : names of patches. Wildcards allowed.
        patches (outlet);
    }

    {
        // Name of new patch
        name wall;

        // Dictionary to construct new patch from
        patchInfo
        {
            type wall;
        }

        // How to construct: either from 'patches' or 'set'
        constructFrom patches;

        // If constructFrom = patches : names of patches. Wildcards allowed.
        patches (wall);
    }
);

// ************************************************************************* //
EOL

verify_file "$createpatch_file"

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

# STEP 4: MESH GENERATION
echo_info "Generating mesh using $mesher mesher..."

if [ "$mesher" == "block" ]; then
    # Use blockMesh for mesh generation
    echo_info "Running blockMesh..."
    blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo_error "blockMesh failed. Check log for details: $case_dir/blockMesh.log"
        cat "$case_dir/blockMesh.log"
        exit 1
    fi
    echo_success "blockMesh completed successfully"
    
    # Make sure mesh is clean and ready
    echo_info "Checking mesh quality..."
    checkMesh -case "$case_dir" -noFunctionObjects > "$case_dir/checkMesh.log" 2>&1
    
    if grep -q "FAILED" "$case_dir/checkMesh.log"; then
        echo_warning "Mesh check found issues, see details in: $case_dir/checkMesh.log"
    else
        echo_success "Mesh quality check passed"
    fi
    
elif [ "$mesher" == "gmsh" ]; then
    # Use Gmsh for mesh generation
    echo_info "Running Gmsh to generate mesh..."
    
    # Make sure Gmsh is available
    if ! command -v gmsh &> /dev/null; then
        echo_error "Gmsh executable not found. Please install Gmsh and add it to your PATH."
        exit 1
    fi
    
    # Generate the mesh using Gmsh
    gmsh -2 "$case_dir/aorta.geo" -o "$case_dir/aorta.msh" -format msh2 -order 1 -v 3 > "$case_dir/gmsh.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo_error "Gmsh failed to generate mesh. Check log: $case_dir/gmsh.log"
        cat "$case_dir/gmsh.log"
        exit 1
    fi
    echo_success "Gmsh mesh generation completed"
    
    # Convert Gmsh mesh to OpenFOAM format
    echo_info "Converting Gmsh mesh to OpenFOAM format..."
    gmshToFoam "$case_dir/aorta.msh" -case "$case_dir" > "$case_dir/gmshToFoam.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo_error "gmshToFoam conversion failed. Check log: $case_dir/gmshToFoam.log"
        cat "$case_dir/gmshToFoam.log"
        exit 1
    fi
    echo_success "Gmsh mesh converted to OpenFOAM format"
    
    # Fix boundary types with createPatch if needed
    echo_info "Updating boundary types..."
    createPatch -overwrite -case "$case_dir" > "$case_dir/createPatch.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo_warning "createPatch had issues, check log: $case_dir/createPatch.log"
    else
        echo_success "Boundary types updated successfully"
    fi
    
    # Check mesh quality
    echo_info "Checking mesh quality..."
    checkMesh -case "$case_dir" -noFunctionObjects > "$case_dir/checkMesh.log" 2>&1
    
    if grep -q "FAILED" "$case_dir/checkMesh.log"; then
        echo_warning "Mesh check found issues, see details in: $case_dir/checkMesh.log"
    else
        echo_success "Mesh quality check passed"
    fi
fi

# STEP 5: INITIALIZE FIELDS ON THE MESH
echo_info "Re-initializing fields on the mesh to ensure consistency..."

# First, check if mesh was generated successfully
if [ ! -d "$case_dir/constant/polyMesh" ] || [ ! -f "$case_dir/constant/polyMesh/points" ]; then
    echo_error "Mesh generation failed. No valid mesh found in $case_dir/constant/polyMesh"
    exit 1
fi

# Copy the 0 directory to 0.orig for backup
mkdir -p "$case_dir/0.orig"
cp "$case_dir/0/"* "$case_dir/0.orig/"

# Use setFields or manual field manipulation to ensure proper initialization
echo_info "Setting initial fields with proper dimensions..."

# Add a verification step to catch mesh-field inconsistencies
for field in p U k omega nut; do
    if [ -f "$case_dir/0/$field" ]; then
        # Run a simple check using foamDictionary to catch basic inconsistencies
        foamDictionary -entry internalField -case "$case_dir" "0/$field" > /dev/null 2>&1
        if [ $? -ne 0 ]; then
            echo_warning "Field $field may have inconsistencies. Reinitializing..."
            # Copy from orig to ensure clean state
            cp "$case_dir/0.orig/$field" "$case_dir/0/$field"
        fi
    fi
done

# STEP 6: RUN THE SIMULATION

# Create a case.foam file for ParaView
touch "$case_dir/case.foam"
echo_success "Created case.foam file for ParaView"

echo_info "Starting OpenFOAM simulation for stenosis level $stenosis_level..."
echo_info "This may take several minutes depending on your computer..."

# Run simpleFoam solver
simpleFoam -case "$case_dir" > "$case_dir/simpleFoam.log" 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo_success "Simulation completed successfully!"
    echo_info "Results available in directory: $case_dir"
    
    # Run post-processing for wall shear stress and sampling
    echo_info "Running post-processing for visualization data..."
    
    # Get the latest time directory
    latest_time=$(find "$case_dir" -maxdepth 1 -name "[0-9]*" | sort -n | tail -1)
    latest_time_value=$(basename "$latest_time")
    
    if [ ! -z "$latest_time_value" ]; then
        echo_info "Using latest time directory: $latest_time_value"
        
        # Calculate wall shear stress
        postProcess -func wallShearStress -case "$case_dir" -time "$latest_time_value" -noFunctionObjects > "$case_dir/postWallShearStress.log" 2>&1
        
        # Sample data along centerline and walls
        postProcess -func sample -case "$case_dir" -time "$latest_time_value" -noFunctionObjects > "$case_dir/postSample.log" 2>&1
        
        # Create basic plots directory
        mkdir -p "plots"
        echo_info "Creating result plots in 'plots' directory..."
        
        if [ -f "plot_results.py" ]; then
            python3 plot_results.py --dirs "$case_dir" --all
            echo_success "Created visualization plots in 'plots' directory"
        else
            echo_warning "plot_results.py not found, skipping automatic plot generation"
        fi
    else
        echo_warning "No time directories found. Post-processing skipped."
    fi
    
    # Create instructions file in the case directory for reference
    cat > "$case_dir/README.txt" << EOL
Aorta Simulation Case: $case_name
Stenosis Level: $stenosis_level
Mesh Method: $mesher
Mesh Quality: $mesh_quality

To visualize this case:
1. Using ParaView: Run 'paraFoam -case $case_dir'
2. Using OpenFOAM utilities: Run 'postProcess -func sample -case $case_dir -latestTime' for sampling data

For extracting wall shear stress data:
Run 'postProcess -func wallShearStress -case $case_dir -latestTime'

For questions or issues, please check the log files in this directory.
EOL
    
    echo_info "You can visualize results with: paraFoam -case $case_dir"
    
else
    echo_error "Simulation encountered errors. Check log file: $case_dir/simpleFoam.log"
    echo_info "Last 20 lines of the log file:"
    tail -n 20 "$case_dir/simpleFoam.log"
    
    # Try to identify specific errors
    if grep -q "attempt to read beyond EOF" "$case_dir/simpleFoam.log"; then
        echo_error "EOF error detected in boundary conditions. This usually indicates a problem with the field files."
    fi
    if grep -q "could not find file" "$case_dir/simpleFoam.log"; then
        echo_error "Missing file error detected. Check that all required files exist."
    fi
    if grep -q "divergence" "$case_dir/simpleFoam.log"; then
        echo_error "Solution divergence detected. Try reducing relaxation factors or using a coarser mesh."
    fi
fi

echo "======================================================="
echo "Simulation complete for stenosis level: $stenosis_level"
echo "======================================================="