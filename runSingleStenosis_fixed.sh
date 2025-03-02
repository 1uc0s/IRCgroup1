#!/bin/bash
# Source OpenFOAM environment
if [ -f "./setup_openfoam_macos.sh" ]; then
    source ./setup_openfoam_macos.sh
else
    echo "❌ ERROR: setup_openfoam_macos.sh not found. Please create this file first."
    exit 1
fi

# Default values
stenosis_level=0.5
case_name="moderate"

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
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--level 0.5] [--name moderate]"
            exit 1
            ;;
    esac
done

# Create directories
case_dir="aorta_simulation"
mkdir -p $case_dir/{0,constant,system}

echo "======================================================="
echo "Processing stenosis level: $case_name ($stenosis_level)"
echo "======================================================="

# First create all OpenFOAM files - MUST BE DONE BEFORE gmshToFoam
echo "Creating OpenFOAM configuration files..."

# Create controlDict
cat > $case_dir/system/controlDict << EOL
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

endTime         2000;  // Iterations for steady-state solver

deltaT          1;     // Time step for steady solver is iteration step

writeControl    timeStep;

writeInterval   200;   // Write results every 200 iterations

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

# Create fvSchemes
cat > $case_dir/system/fvSchemes << EOL
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
cat > $case_dir/system/fvSolution << EOL
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

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
        nSweeps         1;
    }

    "(U|k|omega)Final"
    {
        \$U;
        relTol          0;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;
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
EOL

# Create initial conditions (0 directory)
# p
cat > $case_dir/0/p << EOL
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

# U
cat > $case_dir/0/U << EOL
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

# k
cat > $case_dir/0/k << EOL
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

# omega
cat > $case_dir/0/omega << EOL
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

# Create turbulenceProperties
mkdir -p $case_dir/constant
cat > $case_dir/constant/turbulenceProperties << EOL
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

# Create transportProperties
cat > $case_dir/constant/transportProperties << EOL
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

# Now create Gmsh script with completely corrected physical groups
cat > $case_dir/aorta.geo << EOF
// Parameters for aorta with atherosclerosis
diameter = 25.0;  // Aorta diameter in mm
length = 150.0;   // Length of segment in mm
stenosis_level = $stenosis_level;  // Stenosis level (0.0-1.0)
stenosis_length = 20.0;  // Length of the stenotic region in mm
stenosis_position = length/2;  // Position of stenosis center from inlet
mesh_size = 0.5;  // Default mesh size
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

// Create lines for the 2D mesh
Line(1) = {1, 2};  // Upper wall segment 1
Line(2) = {2, 3};  // Upper wall segment 2
Line(3) = {3, 4};  // Upper wall segment 3
Line(4) = {4, 5};  // Upper wall segment 4
Line(5) = {6, 7};  // Lower wall segment 1
Line(6) = {7, 8};  // Lower wall segment 2
Line(7) = {8, 9};  // Lower wall segment 3
Line(8) = {9, 10}; // Lower wall segment 4
Line(9) = {1, 6};  // Inlet
Line(10) = {5, 10}; // Outlet

// Create the 2D surface
Line Loop(1) = {1, 2, 3, 4, 10, -8, -7, -6, -5, -9};
Plane Surface(1) = {1};

// Extrude to create 3D
out[] = Extrude {0, 0, extrude_depth} {
  Surface{1}; Layers{1}; Recombine;
};

// Define physical groups
// Define boundaries first (for OpenFOAM BC)
Physical Surface("inlet") = {9};
Physical Surface("outlet") = {10};
Physical Surface("wall") = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Surface("frontAndBack") = {1, out[0]};
Physical Volume("fluid") = {out[1]};

// Mesh settings - create structured hex mesh
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8; // Delaunay for quads
Mesh.ElementOrder = 1;
Mesh.Smoothing = 20;
EOF

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
echo "Fixing mesh boundary conditions..."
boundaryFile="$case_dir/constant/polyMesh/boundary"

# Check if boundary file exists
if [ ! -f "$boundaryFile" ]; then
    echo "❌ ERROR: Boundary file not found at $boundaryFile"
    exit 1
fi

# Create a backup of the boundary file
cp "$boundaryFile" "${boundaryFile}.backup"

# Use sed to modify boundary conditions (for macOS compatibility)
sed -i '' -e '/inlet/,/}/s/type.*;/type            patch;/' "$boundaryFile"
sed -i '' -e '/outlet/,/}/s/type.*;/type            patch;/' "$boundaryFile"
sed -i '' -e '/wall/,/}/s/type.*;/type            wall;/' "$boundaryFile"
sed -i '' -e '/frontAndBack/,/}/s/type.*;/type            empty;/' "$boundaryFile"

# Check mesh
echo "Checking mesh quality..."
checkMesh -case $case_dir > $case_dir/checkMesh.log

# Run the simulation
echo "Running OpenFOAM simulation for $case_name stenosis..."
simpleFoam -case $case_dir > $case_dir/simpleFoam.log

echo "✅ Completed simulation for $case_name stenosis"
echo "Results can be visualized using ParaView: paraFoam -case $case_dir"