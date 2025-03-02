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

# Create controlDict
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
        writeControl    writeTime;
        patches         (wall);
    }

    // Calculate y+ to check mesh quality near walls
    yPlus
    {
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeFields     yes;
        writeControl    writeTime;
    }

    // Calculate forces on the wall
    forces
    {
        type            forces;
        libs            (forces);
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
        fields          (".*");
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

# Create initial field files (0 directory)
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
        // inlet settings
    }

    outlet
    {
        // outlet settings 
    }

    wall
    {
        // wall settings
    }

    frontAndBack
    {
        type            empty;
    }
}
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
        // inlet settings
    }

    outlet
    {
        // outlet settings 
    }

    wall
    {
        // wall settings
    }

    frontAndBack
    {
        type            empty;
    }
}
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
        // inlet settings
    }

    outlet
    {
        // outlet settings 
    }

    wall
    {
        // wall settings
    }

    frontAndBack
    {
        type            empty;
    }
}
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
        // inlet settings
    }

    outlet
    {
        // outlet settings 
    }

    wall
    {
        // wall settings
    }

    frontAndBack
    {
        type            empty;
    }
}
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
        // inlet settings
    }

    outlet
    {
        // outlet settings 
    }

    wall
    {
        // wall settings
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

verify_file "$nut_file"

# STEP 2: CREATE GEOMETRY AND MESH

# Create blockMeshDict for the aorta with stenosis
echo_info "Creating blockMesh dictionary..."
block_mesh_file="$case_dir/system/blockMeshDict"
cat > "$block_mesh_file" << EOL
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
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

// Stenosis parameters
radius           12.5;    // Aorta radius in mm
length           150;     // Total length in mm
stenosisLevel    $stenosis_level;   // Level of stenosis (0-1)
stenosisLength   20;      // Length of stenotic region in mm
stenosisPos      75;      // Position of stenosis center

// Calculate narrowed radius at stenosis
stenosisRadius   #calc "radius * (1.0 - stenosisLevel)";

// Number of cells
nCellsLength     100;     // Number of cells along length
nCellsRadial     15;      // Number of cells in radial direction
nCellsThickness  1;       // One cell in z-direction (2D simulation)

// Define z-thickness (small value for quasi-2D)
zThickness       1;

// Define positions
xMin             0;
xInlet           #calc "stenosisPos - stenosisLength/2";
xStenosis        stenosisPos;
xOutlet          #calc "stenosisPos + stenosisLength/2";
xMax             length;

// Proportional cell distributions
firstLength      #calc "(stenosisPos - stenosisLength/2) / length";
stenosisSection  #calc "stenosisLength / length";
lastLength       #calc "1.0 - firstLength - stenosisSection";

vertices
(
    // Bottom face, inlet to stenosis start
    (0        -radius          0)          // 0
    (xInlet   -radius          0)          // 1
    (xInlet   -radius          zThickness) // 2
    (0        -radius          zThickness) // 3
    
    // Bottom face at stenosis
    (xStenosis -stenosisRadius 0)          // 4
    (xStenosis -stenosisRadius zThickness) // 5
    
    // Bottom face, stenosis to outlet
    (xOutlet  -radius          0)          // 6
    (xOutlet  -radius          zThickness) // 7
    (xMax     -radius          0)          // 8
    (xMax     -radius          zThickness) // 9
    
    // Top face, inlet to stenosis start
    (0        radius           0)          // 10
    (xInlet   radius           0)          // 11
    (xInlet   radius           zThickness) // 12
    (0        radius           zThickness) // 13
    
    // Top face at stenosis
    (xStenosis stenosisRadius  0)          // 14
    (xStenosis stenosisRadius  zThickness) // 15
    
    // Top face, stenosis to outlet
    (xOutlet  radius           0)          // 16
    (xOutlet  radius           zThickness) // 17
    (xMax     radius           0)          // 18
    (xMax     radius           zThickness) // 19
);

// Calculate cell distribution
nCellsFirst      #calc "round($firstLength * $nCellsLength)";
nCellsStenosis   #calc "round($stenosisSection * $nCellsLength)";
nCellsLast       #calc "round($lastLength * $nCellsLength)";

// Adjust if rounding created too few/many cells
nCellsCheck      #calc "$nCellsFirst + $nCellsStenosis + $nCellsLast";
nCellsFirst      #calc "$nCellsFirst + ($nCellsLength - $nCellsCheck)";

blocks
(
    // Block 0: Inlet to stenosis start
    hex (0 1 11 10 3 2 12 13)
    ($nCellsFirst $nCellsRadial $nCellsThickness)
    simpleGrading (1 1 1)
    
    // Block 1: Stenosis section
    hex (1 4 14 11 2 5 15 12)
    ($nCellsStenosis $nCellsRadial $nCellsThickness)
    simpleGrading (1 1 1)
    
    // Block 2: Stenosis to outlet
    hex (4 6 16 14 5 7 17 15)
    ($nCellsStenosis $nCellsRadial $nCellsThickness)
    simpleGrading (1 1 1)
    
    // Block 3: Outlet section
    hex (6 8 18 16 7 9 19 17)
    ($nCellsLast $nCellsRadial $nCellsThickness)
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

# Replace the Gmsh call with blockMesh
echo_info "Generating mesh with blockMesh..."
blockMesh -case "$case_dir" > "$case_dir/blockMesh.log" 2>&1

# Exit if blockMesh failed
if [ $? -ne 0 ]; then
    echo_error "blockMesh failed to generate the mesh. Check log for details."
    cat "$case_dir/blockMesh.log"
    exit 1
fi
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

# STEP 3: GENERATE MESH AND RUN SIMULATION
echo_info "Generating 2D mesh with gmsh..."
gmsh -2 "$case_dir/aorta.geo" -o "$case_dir/aorta.msh"

# Exit if gmsh failed
if [ $? -ne 0 ]; then
    echo_error "Gmsh failed to generate the mesh. Exiting."
    exit 1
fi

echo_info "Converting mesh to OpenFOAM format..."
gmshToFoam "$case_dir/aorta.msh" -case "$case_dir"

# Exit if conversion failed
if [ $? -ne 0 ]; then
    echo_error "gmshToFoam failed to convert the mesh. Exiting."
    exit 1
fi

# Fix boundary types
echo_info "Fixing mesh boundary conditions..."
boundaryFile="$case_dir/constant/polyMesh/boundary"

# Check if boundary file exists
if [ ! -f "$boundaryFile" ]; then
    echo_error "Boundary file not found at $boundaryFile"
    echo_info "This may indicate an issue with the mesh generation or conversion."
    echo_info "Check the gmsh and gmshToFoam output for errors."
    ls -la "$case_dir/constant/polyMesh/"
    exit 1
fi

# Create a backup of the boundary file
cp "$boundaryFile" "${boundaryFile}.backup"


# Verify files exist before running
echo_info "Verifying OpenFOAM case setup before running..."
required_files=(
    "$case_dir/0/p"
    "$case_dir/0/U"
    "$case_dir/0/k"
    "$case_dir/0/omega"
    "$case_dir/0/nut"
    "$case_dir/constant/turbulenceProperties"
    "$case_dir/constant/transportProperties"
    "$case_dir/system/controlDict"
    "$case_dir/system/fvSchemes"
    "$case_dir/system/fvSolution"
    "$case_dir/constant/polyMesh/boundary"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo_error "Required file missing: $file"
        exit 1
    else
        echo_info "Verified: $file exists"
    fi
done

# Check mesh
echo_info "Checking mesh quality..."
checkMesh -case "$case_dir" > "$case_dir/checkMesh.log" 2>&1

# Run the simulation
echo_info "Running OpenFOAM simulation for $case_name stenosis..."
simpleFoam -case "$case_dir" > "$case_dir/simpleFoam.log" 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo_success "Completed simulation for $case_name stenosis"
    echo_info "Results can be found in directory: $case_dir"
    echo_info "Visualize results using: paraFoam -case $case_dir"
else
    echo_warning "Simulation may have encountered issues. Check log file: $case_dir/simpleFoam.log"
    # Print the last few lines of the log to help diagnose
    echo_info "Last 10 lines of the log file:"
    tail -n 10 "$case_dir/simpleFoam.log"
fi