#!/bin/bash
# Script to create the OpenFOAM case structure for aorta simulation

# Create base case directory
caseDir="aortaCase"
mkdir -p $caseDir/{0,constant,system}

# Create constant/polyMesh directory (for mesh)
mkdir -p $caseDir/constant/polyMesh

# Create constant/transportProperties for blood properties
mkdir -p $caseDir/constant

echo "Case structure created at: $caseDir"
echo "Next step: Generate mesh with gmsh and convert to OpenFOAM format"