#!/bin/bash
# First source the environment setup script
if [ -f "./setup_openfoam_macos.sh" ]; then
    source ./setup_openfoam_macos.sh
else
    echo "❌ ERROR: setup_openfoam_macos.sh not found. Please create this file first."
    exit 1
fi

# Create base directory
mkdir -p aorta_simulations

# Stenosis levels (as decimal values)
stenosisLevels=(0.0 0.3 0.5 0.7)
stenosisNames=("healthy" "mild" "moderate" "severe")

# Loop through stenosis levels
for i in "${!stenosisLevels[@]}"
do
    level=${stenosisLevels[$i]}
    name=${stenosisNames[$i]}
    
    echo "======================================================="
    echo "Processing stenosis level: $name ($level)"
    echo "======================================================="
    
    # Case directory for this stenosis level
    caseDir="aorta_simulations/stenosis_${name}"
    
    # Create case directory
    mkdir -p $caseDir
    
    # Copy template files
    cp -r aortaCase/constant $caseDir/ 2>/dev/null || mkdir -p $caseDir/constant
    cp -r aortaCase/system $caseDir/ 2>/dev/null || mkdir -p $caseDir/system
    mkdir -p $caseDir/0
    cp -r aortaCase/0/* $caseDir/0/ 2>/dev/null || true
    
    # Create gmsh script with this stenosis level - NOW WITH 3D EXTRUSION
    cat > $caseDir/aorta.geo << EOF
// Parameters for aorta with atherosclerosis
diameter = 25.0;  // Aorta diameter in mm
length = 150.0;   // Length of segment in mm
stenosis_level = $level;  // Stenosis level (0.0-1.0)
stenosis_length = 20.0;  // Length of the stenotic region in mm
stenosis_position = length/2;  // Position of stenosis center from inlet
mesh_size = 0.5;  // Default mesh size
extrude_depth = 1.0;  // Depth for 3D extrusion (thin layer for 2D simulation)

// Calculate stenosis height
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

// Create splines for upper and lower walls
Spline(1) = {1, 2, 3, 4, 5};  // Upper wall
Spline(2) = {6, 7, 8, 9, 10};  // Lower wall

// Create inlet and outlet lines
Line(3) = {1, 6};  // Inlet
Line(4) = {5, 10};  // Outlet

// Create surface
Line Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};

// Define front and back faces for 3D extrusion
// Extrude the 2D surface in the Z direction to create a 3D mesh
mesh_tags[] = Extrude(0, 0, extrude_depth) {
  Surface{1}; Layers{1}; Recombine;
};

// Physical groups for boundary conditions
Physical Surface("inlet") = {3};  // Inlet face becomes a surface in 3D
Physical Surface("outlet") = {4}; // Outlet face becomes a surface in 3D
Physical Surface("wall") = {1, 2}; // Wall faces become surfaces in 3D
Physical Surface("frontAndBack") = {1, mesh_tags[0]}; // Front and back faces
Physical Volume("fluid") = {mesh_tags[1]}; // 3D volume

// Mesh settings
Mesh.RecombineAll = 1;  // Generate quadrilateral/hexahedral elements
Mesh.Smoothing = 20;    // Mesh smoothing steps
EOF
    
    echo "Generating 3D mesh with gmsh for $name stenosis..."
    gmsh -3 $caseDir/aorta.geo -o $caseDir/aorta.msh
    
    # Exit if gmsh failed
    if [ $? -ne 0 ]; then
        echo "❌ ERROR: Gmsh failed to generate the mesh. Exiting."
        exit 1
    fi
    
    echo "Converting mesh to OpenFOAM format..."
    gmshToFoam $caseDir/aorta.msh -case $caseDir
    
    # Exit if conversion failed
    if [ $? -ne 0 ]; then
        echo "❌ ERROR: gmshToFoam failed to convert the mesh. Exiting."
        exit 1
    fi
    
    # Fix boundary types
    echo "Fixing mesh boundary conditions..."
    boundaryFile="$caseDir/constant/polyMesh/boundary"
    
    # Check if boundary file exists
    if [ ! -f "$boundaryFile" ]; then
        echo "❌ ERROR: Boundary file not found at $boundaryFile"
        exit 1
    fi
    
    # Using perl for macOS compatibility
    perl -i -pe 's/(inlet.*?type\s+)(\w+)(;)/\1patch\3/gs' "$boundaryFile"
    perl -i -pe 's/(outlet.*?type\s+)(\w+)(;)/\1patch\3/gs' "$boundaryFile"
    perl -i -pe 's/(wall.*?type\s+)(\w+)(;)/\1wall\3/gs' "$boundaryFile"
    perl -i -pe 's/(frontAndBack.*?type\s+)(\w+)(;)/\1empty\3/gs' "$boundaryFile"
    
    # Create the basic OpenFOAM configuration files if they don't exist
    echo "Creating OpenFOAM configuration files..."
    
    # Create controlDict if it doesn't exist
    if [ ! -f "$caseDir/system/controlDict" ]; then
        mkdir -p $caseDir/system
        cat > $caseDir/system/controlDict << EOL
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

endTime         5000;  // Iterations for steady-state solver

deltaT          1;     // Time step for steady solver is iteration step

writeControl    timeStep;

writeInterval   500;   // Write results every 500 iterations

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
}
EOL
    fi
    
    # Create fvSchemes if it doesn't exist
    if [ ! -f "$caseDir/system/fvSchemes" ]; then
        mkdir -p $caseDir/system
        cat > $caseDir/system/fvSchemes << EOL
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
    fi
    
    # Create fvSolution if it doesn't exist
    if [ ! -f "$caseDir/system/fvSolution" ]; then
        mkdir -p $caseDir/system
        cat > $caseDir/system/fvSolution << EOL
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
    fi
    
    # Create initial conditions (0 directory)
    # p
    if [ ! -f "$caseDir/0/p" ]; then
        mkdir -p $caseDir/0
        cat > $caseDir/0/p << EOL
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
    fi
    
    # U
    if [ ! -f "$caseDir/0/U" ]; then
        mkdir -p $caseDir/0
        cat > $caseDir/0/U << EOL
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
    fi
    
    # k
    if [ ! -f "$caseDir/0/k" ]; then
        mkdir -p $caseDir/0
        cat > $caseDir/0/k << EOL
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
    fi
    
    # omega
    if [ ! -f "$caseDir/0/omega" ]; then
        mkdir -p $caseDir/0
        cat > $caseDir/0/omega << EOL
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
    fi
    
    # Create turbulenceProperties
    if [ ! -f "$caseDir/constant/turbulenceProperties" ]; then
        mkdir -p $caseDir/constant
        cat > $caseDir/constant/turbulenceProperties << EOL
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
    fi
    
    # Create transportProperties
    if [ ! -f "$caseDir/constant/transportProperties" ]; then
        mkdir -p $caseDir/constant
        cat > $caseDir/constant/transportProperties << EOL
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
    fi
    
    # Check mesh
    echo "Checking mesh quality..."
    checkMesh -case $caseDir > $caseDir/checkMesh.log
    
    # Run the simulation
    echo "Running OpenFOAM simulation for $name stenosis..."
    simpleFoam -case $caseDir > $caseDir/simpleFoam.log
    
    echo "✅ Completed simulation for $name stenosis"
    echo ""
done

echo "🎉 All simulations completed!"
echo "Results can be visualized using ParaView."