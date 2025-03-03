#!/bin/bash
# runStenosis_highStenosis.sh - Specialized script for high stenosis cases
# Uses a two-stage approach: steady-state initialization followed by pulsatile flow

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
case_name="severe_pulsatile"
mesh_size=1.0   # Coarser mesh for initial stability
time_step=0.001
num_cycles=5
cycle_time=1.0
init_time=2.0   # Time for steady-state initialization

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
        --time-step)
            time_step="$2"
            shift 2
            ;;
        --cycles)
            num_cycles="$2"
            shift 2
            ;;
        --cycle-time)
            cycle_time="$2"
            shift 2
            ;;
        --init-time)
            init_time="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --level VALUE      Stenosis level (0.0-0.9), default: 0.5"
            echo "  --name NAME        Name for the case, default: severe_pulsatile"
            echo "  --mesh SIZE        Mesh size (larger = coarser), default: 1.0"
            echo "  --time-step TIME   Time step size in seconds, default: 0.001"
            echo "  --cycles NUM       Number of cardiac cycles, default: 5"
            echo "  --cycle-time TIME  Duration of one cardiac cycle in seconds, default: 1.0"
            echo "  --init-time TIME   Time for steady-state initialization, default: 2.0"
            echo "  --help             Display this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run '$0 --help' for usage information."
            exit 1
            ;;
    esac
done

# Create case directory
case_dir="aorta_simulation_${case_name}"
echo_info "Creating simulation directory structure: $case_dir"

# Create directories
mkdir -p "$case_dir"
mkdir -p "$case_dir/0"
mkdir -p "$case_dir/0.orig"
mkdir -p "$case_dir/constant"
mkdir -p "$case_dir/system"

# Verify directories were created
for dir in "$case_dir" "$case_dir/0" "$case_dir/0.orig" "$case_dir/constant" "$case_dir/system"; do
    if [ -d "$dir" ]; then
        echo_success "Created directory: $dir"
    else
        echo_error "Failed to create directory: $dir"
        exit 1
    fi
done

echo "======================================================="
echo "Processing high stenosis (${stenosis_level}) with pulsatile flow"
echo "Using two-stage approach: steady-state initialization followed by pulsatile"
echo "======================================================="

# Calculate the total simulation time
total_time=$(echo "$init_time + ($num_cycles * $cycle_time)" | bc -l)

# Create controlDict for the two-stage simulation
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

// Initially using PISO with simple settings for stability
application     pisoFoam;   

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         $total_time;

deltaT          $time_step;

writeControl    runTime;

writeInterval   0.1;   // Write every 0.1s

purgeWrite      0;     // Keep all time directories

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

// Define a solver switch at initialization point
functions
{
    // Switch from PISO to PIMPLE at a specific time
    solverSwitch
    {
        type               coded;
        name               switchSolver;
        writeControl       timeStep;
        writeInterval      1;
        
        codeExecute
        #{
            // Get current time
            scalar currentTime = this->mesh().time().value();
            
            if (currentTime >= $init_time && currentTime <= ($init_time + 0.1))
            {
                Info<< "====================================================" << endl;
                Info<< "Switching from initialization to pulsatile flow at time " << currentTime << endl;
                Info<< "====================================================" << endl;
                
                // No actual action needed here - just a marker for post-processing
            }
        #};
    }
    
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

# Create fvSchemes with first-order schemes for stability
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

// First-order schemes for stability
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default          none;
    div(phi,U)       Gauss upwind;  // First-order upwind for stability
    div(phi,k)       Gauss upwind;
    div(phi,omega)   Gauss upwind;
    div(phi,nut)     Gauss upwind;
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

# Create fvSolution with PISO algorithm settings - simplest possible for now
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
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-6;
        relTol          0.1;
    }
    
    pFinal
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-6;
        relTol          0.1;
    }

    "(U|k|omega|nut)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }

    "(U|k|omega|nut)Final"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }
}

PISO
{
    nCorrectors      2;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.5;
        k               0.5;
        omega           0.5;
    }
}
EOL

verify_file "$fv_solution_file"

# Create turbulenceProperties with SST model
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

# Create initial conditions - simple steady inflow
# p (pressure)
p_file="$case_dir/0.orig/p"
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

verify_file "$p_file"

# Copy to 0 directory
cp "$p_file" "$case_dir/0/"

# U (velocity) - Low steady velocity for initialization
u_file="$case_dir/0.orig/U"
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

internalField   uniform (0.05 0 0);  // Very low velocity for stability

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (0.05 0 0);
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

verify_file "$u_file"

# Copy to 0 directory
cp "$u_file" "$case_dir/0/"

# Create a special pulsatile velocity file for later use
pulsatile_u_file="$case_dir/0.orig/U.pulsatile"
cat > "$pulsatile_u_file" << EOL
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

internalField   uniform (0.05 0 0);  // Keep velocity low

boundaryField
{
    inlet
    {
        type            codedFixedValue;
        value           uniform (0.05 0 0);
        
        name            pulsatileFlow;
        
        code
        #{
            // Calculate phase in cardiac cycle (0 to 1)
            scalar t = this->db().time().value();
            scalar initTime = $init_time;
            scalar cycleTime = $cycle_time;
            
            // Only apply pulsatile flow after initialization period
            if (t <= initTime) {
                // During initialization, keep constant low velocity
                scalar velocity = 0.05;
                
                const vectorField& centers = patch().Cf();
                vectorField& patchField = *this;
                
                forAll(patchField, i) {
                    patchField[i] = vector(velocity, 0, 0);
                }
                
                return;
            }
            
            // Calculate phase for pulsatile period
            scalar adjustedTime = t - initTime;
            scalar phase = fmod(adjustedTime, cycleTime) / cycleTime;
            
            // Use very low velocities for stability with high stenosis
            scalar baseVelocity = 0.05;  // m/s, minimum velocity
            scalar amplitude = 0.05;     // m/s, very small amplitude for stability
            
            // Calculate pulsatile velocity using a very smooth profile
            scalar velocity;
            
            // Create a very smooth transition for better numerical stability
            if (phase < 0.3) {
                // Systolic acceleration (smoother cosine rise)
                scalar normPhase = phase / 0.3;  // normalize to 0-1 for this segment
                velocity = baseVelocity + amplitude * 0.5 * (1.0 - cos(normPhase * M_PI));
            } else if (phase < 0.6) {
                // Systolic deceleration (smoother fall)
                scalar normPhase = (phase - 0.3) / 0.3;  // normalize to 0-1 for this segment
                velocity = baseVelocity + amplitude * 0.5 * (1.0 + cos(normPhase * M_PI));
            } else {
                // Diastolic phase (low steady flow)
                velocity = baseVelocity;
            }
            
            // Get mesh points for inlet patch
            const vectorField& centers = patch().Cf();
            
            // Apply parabolic profile across the inlet (assuming circular inlet)
            const scalar radius = 0.0125;  // 12.5mm radius of aorta in meters
            
            vectorField& patchField = *this;
            
            forAll(centers, i) {
                vector pos = centers[i];
                
                // Calculate distance from center of inlet
                scalar y = pos.y();
                scalar z = pos.z();
                scalar distFromCenter = sqrt(y*y + z*z);
                
                // Apply parabolic profile: v = v_max * (1 - (r/R)^2)
                scalar velFactor = 1.0 - pow(distFromCenter/radius, 2);
                if (velFactor < 0) velFactor = 0;
                
                patchField[i] = vector(velocity * velFactor, 0, 0);
            }
        #};
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

verify_file "$pulsatile_u_file"

# Create k (turbulent kinetic energy)
k_file="$case_dir/0.orig/k"
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

internalField   uniform 0.0001;  // Lower initial value for stability

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0.0001;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            kqRWallFunction;
        value           uniform 0.0001;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

verify_file "$k_file"

# Copy to 0 directory
cp "$k_file" "$case_dir/0/"

# Create omega (specific dissipation rate)
omega_file="$case_dir/0.orig/omega"
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

internalField   uniform 10;  // Lower value for stability

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 10;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            omegaWallFunction;
        value           uniform 10;
    }

    frontAndBack
    {
        type            empty;
    }
}
EOL

verify_file "$omega_file"

# Copy to 0 directory
cp "$omega_file" "$case_dir/0/"

# Create nut (turbulent viscosity)
nut_file="$case_dir/0.orig/nut"
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

verify_file "$nut_file"

# Copy to 0 directory
cp "$nut_file" "$case_dir/0/"

# Create blockMeshDict
echo_info "Creating blockMesh dictionary..."
block_mesh_file="$case_dir/system/blockMeshDict"

# Calculate stenosis radius
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

# Calculate cell counts - use coarser mesh for stability
TOTAL_CELLS_LENGTH=80  # Reduced from 100
CELLS_FIRST=$(echo "($FIRST_LENGTH * $TOTAL_CELLS_LENGTH)/1" | bc)
CELLS_STENOSIS=$(echo "($STENOSIS_SECTION * $TOTAL_CELLS_LENGTH)/1" | bc)
CELLS_LAST=$(echo "($LAST_LENGTH * $TOTAL_CELLS_LENGTH)/1" | bc)
CELLS_CHECK=$(echo "$CELLS_FIRST + $CELLS_STENOSIS + $CELLS_LAST" | bc)

# Adjust if rounding created too few/many cells
if [ "$CELLS_CHECK" -ne "$TOTAL_CELLS_LENGTH" ]; then
    CELLS_FIRST=$(echo "$CELLS_FIRST + ($TOTAL_CELLS_LENGTH - $CELLS_CHECK)" | bc)
fi

# Create blockMeshDict with explicitly calculated values
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

// Stenosis parameters
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
    ($CELLS_FIRST 12 1)  // Coarser mesh
    simpleGrading (1 1 1)
    
    // Block 1: Stenosis section (first half)
    hex (1 4 14 11 2 5 15 12)
    ($CELLS_STENOSIS 12 1)  // Coarser mesh
    simpleGrading (1 1 1)
    
    // Block 2: Stenosis section (second half)
    hex (4 6 16 14 5 7 17 15)
    ($CELLS_STENOSIS 12 1)  // Coarser mesh
    simpleGrading (1 1 1)
    
    // Block 3: Outlet section
    hex (6 8 18 16 7 9 19 17)
    ($CELLS_LAST 12 1)  // Coarser mesh
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

# Generate mesh
echo_info "Generating mesh..."
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

# Run the simulation in multiple stages

# Stage 1: Run initialization with PISO (steady-state mode)
echo_info "Starting initialization phase with pisoFoam for $init_time seconds..."
pisoFoam -case "$case_dir" > "$case_dir/pisoFoam.log" 2>&1 &
PISO_PID=$!

# Monitor the progress of the initialization
echo_info "Monitoring initialization progress..."
PREV_TIME=0
while kill -0 $PISO_PID 2>/dev/null; do
    if [ -d "$case_dir" ]; then
        LATEST_TIME=$(find "$case_dir" -maxdepth 1 -name "[0-9]*" | sort -n | tail -1)
        if [ -n "$LATEST_TIME" ]; then
            CURRENT_TIME=$(basename "$LATEST_TIME")
            if (( $(echo "$CURRENT_TIME > $PREV_TIME" | bc -l) )); then
                echo_info "Simulation time: $CURRENT_TIME seconds"
                PREV_TIME=$CURRENT_TIME
                
                # Check if we've reached initialization time
                if (( $(echo "$CURRENT_TIME >= $init_time" | bc -l) )); then
                    echo_info "Initialization phase complete, switching to pulsatile flow..."
                    break
                fi
            fi
        fi
    fi
    sleep 5
done

# If pisoFoam is still running, stop it
if kill -0 $PISO_PID 2>/dev/null; then
    echo_info "Stopping pisoFoam to switch to pulsatile flow..."
    kill $PISO_PID
    wait $PISO_PID 2>/dev/null
fi

# Stage 2: Switch to pulsatile flow
echo_info "Applying pulsatile flow boundary conditions..."

# Find the latest time directory
LATEST_TIME=$(find "$case_dir" -maxdepth 1 -name "[0-9]*" | sort -n | tail -1)
if [ -z "$LATEST_TIME" ]; then
    echo_error "No time directories found. Initialization phase failed."
    exit 1
fi

CURRENT_TIME=$(basename "$LATEST_TIME")
echo_info "Using time directory: $CURRENT_TIME"

# Switch to pulsatile velocity boundary condition
cp "$case_dir/0.orig/U.pulsatile" "$LATEST_TIME/U"
echo_success "Applied pulsatile boundary condition at time $CURRENT_TIME"

# Change the application in controlDict to pimpleFoam
sed -i '' 's/application     pisoFoam/application     pimpleFoam/' "$case_dir/system/controlDict"

# Update fvSolution to use PIMPLE algorithm
cat > "$case_dir/system/fvSolution" << EOL
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
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-6;
        relTol          0.1;
    }
    
    pFinal
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-6;
        relTol          0.1;
    }

    "(U|k|omega|nut)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }

    "(U|k|omega|nut)Final"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }
}

PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors      1;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
    momentumPredictor no;
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.3;
        k               0.3;
        omega           0.3;
    }
}
EOL

# Run the pulsatile simulation
echo_info "Starting pulsatile flow simulation with pimpleFoam..."
echo_info "Simulating from time $CURRENT_TIME to $total_time seconds..."

pimpleFoam -case "$case_dir" > "$case_dir/pimpleFoam.log" 2>&1

# Check if simulation completed successfully
if [ $? -eq 0 ]; then
    echo_success "Simulation completed successfully!"
    echo_info "Results available in directory: $case_dir"
    
    # Create basic plots directory
    mkdir -p "plots"
    echo_info "Creating basic result plots in 'plots' directory..."
    
    if [ -f "plot_pulsatile_results.py" ]; then
        python3 plot_pulsatile_results.py --case "$case_dir" --all
        echo_success "Created visualization plots in 'plots' directory"
    elif [ -f "plot_results.py" ]; then
        python3 plot_results.py --dirs "$case_dir" --all
        echo_success "Created basic visualization plots in 'plots' directory"
    else
        echo_warning "Plotting scripts not found, skipping automatic plot generation"
    fi
    
    echo_info "For animation of pulsatile flow, use: paraFoam -case $case_dir"
else
    echo_error "Simulation encountered errors. Check log file: $case_dir/pimpleFoam.log"
    echo_info "Last 20 lines of the log file:"
    tail -n 20 "$case_dir/pimpleFoam.log"
fi

echo "======================================================="
echo "Two-stage pulsatile flow simulation complete for stenosis level: $stenosis_level"
echo "Initialization time: $init_time seconds"
echo "Pulsatile cycles: $num_cycles with cycle time of $cycle_time seconds"
echo "======================================================="