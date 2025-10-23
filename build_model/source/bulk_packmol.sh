#!/bin/bash

# Script to automate system setup using Packmol and GROMACS
set -euo pipefail
echo -e "${BLUE}
  [PIPELINE] Starting bulk preparation
  Script    : ${BASH_SOURCE[0]}
${NC}" | tee -a ./result/b_model.log >&2
cp "${bulk_pa}".inp bulk_.inp
# Step 1: Run Packmol to generate initial configuration

echo "Running Packmol to generate PDB file..." | tee -a ./result/b_model.log >&2
if packmol < bulk_.inp >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}Packmol completed successfully." | tee -a ./result/b_model.log >&2
else
    echo -e "${RED}Error: Packmol failed. Check ${bulk_pa}.inp for errors.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

bulk_name=$(sed -n 's/^output \([^.]*\)\.pdb.*/\1/p' bulk_.inp)
# Step 3: Convert PDB to GROMACS GRO format using editconf
if gmx editconf -f "${bulk_name}".pdb -o "${bulk_name}s".gro >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}Structure converted to GRO format." | tee -a ./result/b_model.log >&2
else
    echo -e "${RED}Error: gmx editconf failed. Is GROMACS installed and in PATH?${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

echo q | gmx make_ndx -f "${bulk_name}s".gro >> ./result/b_model.log 2>&1

echo "Setting box vectors in "${bulk_name}s".gro to: ${bulk_xbox} x ${bulk_ybox} x ${bulk_zbox} nm" | tee -a ./result/b_model.log >&2
formatted_line=$(printf "%8.4f%8.4f%8.4f" $bulk_xbox $bulk_ybox $bulk_zbox)
sed -i "\$ s/^.*$/      $formatted_line/" "${bulk_name}s".gro
echo -e "${OK}${GREEN}Pipeline completed successfully!${NC}" | tee -a ./result/b_model.log >&2
echo "Final structure: "${bulk_name}s".gro"

rm "${bulk_name}".pdb