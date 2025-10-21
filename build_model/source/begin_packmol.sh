#!/bin/bash
# Script to automate system setup using Packmol and GROMACS
set -euo pipefail
####### add 
echo -e "${BLUE}
 ---------------------------------------- iter="${iter}" ----------------------------------------- 
  [PIPELINE] Starting system preparation
  Script    : ${BASH_SOURCE[0]}
${NC}" | tee -a ./result/b_model.log >&2


# Step 1: Run Packmol to generate initial configuration
echo "Running Packmol to generate PDB file..." | tee -a ./result/b_model.log >&2
if packmol < packmol_.inp >> ./result/b_model.log 2>&1; then
    echo -e "${OK}Packmol completed successfully." | tee -a ./result/b_model.log >&2
else
    echo -e "${RED}Error: Packmol failed. Check packmol.inp for errors.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# Step 2: Extract output filename from packmol.inp
name=$(sed -n 's/^output \([^.]*\)\.pdb.*/\1/p' packmol_.inp)
if [ -z "$name" ]; then
    echo -e "${RED}Error: Could not extract output name from packmol.inp. Make sure 'output *.pdb' is defined.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# Step 3: Convert PDB to GROMACS GRO format using editconf
if gmx editconf -f "$name".pdb -o "${name}s".gro >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}Structure converted to GRO format.">&2
else
    echo -e "${RED}Error: gmx editconf failed. Is GROMACS installed and in PATH?${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# Step 4: Clean up backup files
rm -f *#

# Step 5: Create index file (optional step; here we just run make_ndx and quit)
echo q | gmx make_ndx -f "${name}s".gro >> ./result/b_model.log 2>&1

# Step 6: (Optional) Update box vectors in the .gro file
echo "Setting box vectors in ${name}s.gro to: $xbox x $ybox x $zbox nm" | tee -a ./result/b_model.log >&2
formatted_line=$(printf "%8.4f%8.4f%8.4f" $xbox $ybox $zbox)
sed -i "\$ s/^.*$/      $formatted_line/" "${name}s".gro
echo -e "${OK}${GREEN}Pipeline completed successfully!${NC}" | tee -a ./result/b_model.log >&2
echo "Final structure: ${name}s.gro"

rm ${name}.pdb