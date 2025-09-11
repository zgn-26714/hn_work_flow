#!/bin/bash

# Script to automate system setup using Packmol and GROMACS
set -euo pipefail
echo "
===================================================================================================
 ---------------------------------------- iter="${iter}" -----------------------------------------
  [PIPELINE] Starting system preparation
  Script    : ${BASH_SOURCE[0]}
===================================================================================================">&2

# Step 1: Run Packmol to generate initial configuration
echo "Running Packmol to generate PDB file...">&2
if packmol < packmol.inp > test.log 2>&1 ; then
    echo "Packmol completed successfully.">&2
else
    echo "Error: Packmol failed. Check packmol.inp for errors." >&2
    exit 1
fi

# Step 2: Extract output filename from packmol.inp
name=$(sed -n 's/^output \([^.]*\)\.pdb.*/\1/p' packmol.inp)
if [ -z "$name" ]; then
    echo "Error: Could not extract output name from packmol.inp. Make sure 'output *.pdb' is defined." >&2
    exit 1
fi

# Step 3: Convert PDB to GROMACS GRO format using editconf
if gmx editconf -f "$name".pdb -o "$name".gro > test.log 2>&1 ; then
    echo "Structure converted to GRO format.">&2
else
    echo "Error: gmx editconf failed. Is GROMACS installed and in PATH?" >&2
    exit 1
fi

# Step 4: Clean up backup files
rm -f *#

# Step 5: Create index file (optional step; here we just run make_ndx and quit)
echo q | gmx make_ndx -f "$name".gro > test.log 2>&1

# Step 6: (Optional) Update box vectors in the .gro file
echo "Setting box vectors in ${name}.gro to: $xbox x $ybox x $zbox nm">&2
sed -i '$ s/^.*$/        "$xbox".0000   "$ybox".0000   "$zbox".0000             /' "$name".gro
echo "Pipeline completed successfully!">&2
echo "Final structure: ${name}.gro"