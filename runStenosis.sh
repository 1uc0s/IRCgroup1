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

# Create Gmsh script with parametric stenosis
echo_info "Creating Gmsh geometry script..."
geo_file="$case_dir/aorta.geo"
cat > "$geo_file" << EOF
// Parameters for aorta with atherosclerosis
diameter = 25.0;  // Aorta diameter in mm
length = 150.0;   // Length of segment in mm
stenosis_level = $stenosis_level;  // Stenosis level (0.0-1.0)
stenosis_length = 20.0;  // Length of the stenotic region in mm
stenosis_position = length/2;  // Position of stenosis center from inlet
mesh_size = $mesh_size;  // Default mesh size
extrude_depth = 1.0;  // Depth for 3D extrusion

// Calculate stenosis height based on stenosis level
stenosis_height = diameter * stenosis_level / 2;

// Points for upper wall
Point(1) = {0, diameter/2, 0, mesh_size};  // Inlet top
Point(2) = {stenosis_position - stenosis_length/2, diameter/2, 0, mesh_size};  // Start of stenosis top
Point(3) = {stenosis_position, diameter/2 - stenosis_height, 0, mesh_size/4};  // Peak of stenosis top 
Point(4) = {stenosis_position + stenosis_length/2, diameter/2, 0, mesh_size};  // End of stenosis top
Point(5) = {length, diameter/2, 0, mesh_size};  // Outlet top

// Points for lower wall
Point(6) = {0, -diameter/2, 0, mesh_size};  // Inlet bottom
Point(7) = {stenosis_position - stenosis_length/2, -diameter/2, 0, mesh_size};  // Start of stenosis bottom
Point(8) = {stenosis_position, -diameter/2 + stenosis_height, 0, mesh_size/4};  // Peak of stenosis bottom 
Point(9) = {stenosis_position + stenosis_length/2, -diameter/2, 0, mesh_size};  // End of stenosis bottom
Point(10) = {length, -diameter/2, 0, mesh_size};  // Outlet bottom

// Create splines for the walls (smoother than straight lines)
Spline(1) = {1, 2, 3, 4, 5};  // Upper wall
Spline(2) = {6, 7, 8, 9, 10};  // Lower wall

// Create inlet and outlet lines
Line(3) = {1, 6};  // Inlet
Line(4) = {5, 10}; // Outlet

// Create the 2D surface
Line Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};

// Extrude to create 3D with single layer
// This creates a thin 3D mesh with quad elements
out[] = Extrude {0, 0, extrude_depth} {
  Surface{1}; Layers{1}; Recombine;
};

// Define physical entities with clear naming
Physical Surface("inlet") = {3};       // Inlet face
Physical Surface("outlet") = {4};      // Outlet face
Physical Surface("wall") = {1, 2};     // Wall faces (upper and lower)
Physical Surface("frontAndBack") = {1, out[0]}; // Front and back faces
Physical Volume("fluid") = {out[1]};   // Fluid volume

// Mesh settings - create structured hex mesh
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8; // Delaunay for quads
Mesh.ElementOrder = 1;
Mesh.Smoothing = 20;
EOF

verify_file "$geo_file"

echo "Generating 3D mesh with gmsh..."
gmsh -3 $case_dir/aorta.geo -o $case_dir/aorta.msh

# Exit if gmsh failed
if [ $? -ne 0 ]; then
    echo "❌ ERROR: Gmsh failed to generate the mesh. Exiting."
    exit 1
fi

echo "Converting mesh to OpenFOAM format..."
gmshToFoam $case_dir/aorta.msh -case $case_dir

# Exit if conversion failed
if [ $? -ne 0 ]; then
    echo "❌ ERROR: gmshToFoam failed to convert the mesh. Exiting."
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

# Use sed to modify boundary conditions (compatible with macOS and Linux)
if [ "$OS" == "Darwin" ]; then
    # macOS version of sed requires empty string after -i
    sed -i '' -e '/inlet/,/}/s/type.*;/type            patch;/' "$boundaryFile"
    sed -i '' -e '/outlet/,/}/s/type.*;/type            patch;/' "$boundaryFile"
    sed -i '' -e '/wall/,/}/s/type.*;/type            wall;/' "$boundaryFile"
    sed -i '' -e '/frontAndBack/,/}/s/type.*;/type            empty;/' "$boundaryFile"
else
    # Linux version of sed
    sed -i -e '/inlet/,/}/s/type.*;/type            patch;/' "$boundaryFile"
    sed -i -e '/outlet/,/}/s/type.*;/type            patch;/' "$boundaryFile"
    sed -i -e '/wall/,/}/s/type.*;/type            wall;/' "$boundaryFile"
    sed -i -e '/frontAndBack/,/}/s/type.*;/type            empty;/' "$boundaryFile"
fi

# Create OpenFOAM configuration files
echo "Creating OpenFOAM configuration files..."

# Create controlDict
echo_info "Creating OpenFOAM configuration files..."
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

verify_file "$nut_file"

verify_file "$omega_file"

verify_file "$k_file"

verify_file "$u_file"

verify_file "$p_file"

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

# Create initial conditions (0 directory)
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
EOL

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
EOL

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
        value           \$internalField;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            kqRWallFunction;
        value           \$internalField;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

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
        value           \$internalField;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            omegaWallFunction;
        value           \$internalField;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# nut (turbulent viscosity) - THIS WAS MISSING BEFORE
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
        value           \$internalField;
    }

    outlet
    {
        type            calculated;
        value           \$internalField;
    }

    wall
    {
        type            nutkWallFunction;
        value           \$internalField;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

# Create turbulenceProperties
echo_info "Creating physical properties files..."
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

# Check mesh
echo "Checking mesh quality..."
checkMesh -case $case_dir > $case_dir/checkMesh.log

# Verify all required files exist before running the simulation
echo_info "Verifying all required files exist..."
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
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo_error "Required file missing: $file"
        echo_info "Trying to create this file now..."
        
        # Determine which file is missing and recreate it
        case "$file" in
            *p)
                cat > "$file" << EOL
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
                ;;
            *U)
                # Similar blocks for each required file
                # ... (other file recreation logic)
                ;;
        esac
        
        if [ -f "$file" ]; then
            echo_success "Created missing file: $file"
        else
            echo_error "Failed to create file: $file"
            exit 1
        fi
    fi
done

# Run the simulation
echo_info "Running OpenFOAM simulation for $case_name stenosis..."
simpleFoam -case $case_dir > $case_dir/simpleFoam.log 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo "✅ Completed simulation for $case_name stenosis"
    echo "Results can be found in directory: $case_dir"
    echo "Visualize results using: paraFoam -case $case_dir"
else
    echo "⚠️ Simulation may have encountered issues. Check log file: $case_dir/simpleFoam.log"
fi