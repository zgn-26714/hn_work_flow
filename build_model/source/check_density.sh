#!/bin/bash

# =============================
# Safety settings and checks
# =============================

set -euo pipefail
echo -e "${BLUE}
 ---------------------------------------- iter="${iter}" -----------------------------------------
  "[INFO]Starting pre-equilibration and density analysis workflow"
  Script    : ${BASH_SOURCE[0]}
${NC}"  | tee -a ./result/b_model.log >&2


# Clean up old files safely
echo "Cleaning up temporary files..." | tee -a ./result/b_model.log >&2
rm -f *#
mkdir -p build
echo "[STEP 1]Running grompp to generate mini_pre_eq.tpr" | tee -a ./result/b_model.log >&2

# make index file
echo q | gmx make_ndx -f "${name_value}.gro" -o ./build/index.ndx >> ./result/b_model.log 2>&1
echo -e "${OK}success generate index." | tee -a ./result/b_model.log >&2

if gmx grompp \
    -f "${mini_MDP}.mdp" \
    -c "${name_value}.gro" \
    -p "${TOP}.top" \
    -n ./build/index.ndx \
    -o ./build/mini_pre_eq.tpr \
    -maxwarn "${maxWarn}" >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}grompp completed successfully" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}grompp failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# =============================
# Step 2: mdrun - run simulation
# =============================
echo "[STEP 2]Running mdrun for pre-equilibration" | tee -a ./result/b_model.log >&2
if gmx mdrun \
    -s ./build/mini_pre_eq.tpr \
    -deffnm ./build/mini_pre_eq \
    -ntmpi 1 \
    -ntomp "$NPOS" \
    -v >> ./result/b_model.log 2>&1 ; then
    echo -e "${OK}min energy completed successfully" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}min energy failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# Rename output file
mv ./build/mini_pre_eq.gro ./build/pre_eq.gro


# =============================
# Step 1: grompp - generate tpr
# =============================
echo "[STEP 3]Running grompp to generate pre_eq.tpr" | tee -a ./result/b_model.log >&2

cp "${MDP}".mdp ./build/grompp.mdp
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_den}/" ./build/grompp.mdp
if gmx grompp \
    -f ./build/grompp.mdp \
    -c ./build/pre_eq.gro \
    -p "${TOP}.top" \
    -n ./build/index.ndx \
    -o ./build/pre_eq.tpr \
    -maxwarn "${maxWarn}" >> ./result/b_model.log  2>&1; then
    echo -e "${OK}grompp completed successfully" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}grompp failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# =============================
# Step 2: mdrun - run simulation
# =============================
echo "[STEP 4]Running mdrun for pre-equilibration" | tee -a ./result/b_model.log >&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\tUsing GPU acceleration" | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./build/pre_eq.tpr
        -deffnm ./build/pre_eq
        -ntmpi 1
        -ntomp "$NPOS"
        -pme gpu
        -pmefft gpu
        -nb gpu
        -tunepme no
        -v
    )
else
    echo -e "\tUsing CPU computation" | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./build/pre_eq.tpr
        -deffnm ./build/pre_eq
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi
"${mdrun_cmd[@]}" >> ./result/b_model.log  2>&1 

# Execute mdrun
if [[ -s ./build/pre_eq.xtc ]]; then
    echo -e "${OK}mdrun completed successfully" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}mdrun failed: pre_eq.edr not generated or empty${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# =============================
# Step 3: Compute density profile
# =============================
echo "[STEP 5]Computing density profile" | tee -a ./result/b_model.log >&2
if echo 0|gmx density \
    -f ./build/pre_eq.xtc \
    -s ./build/pre_eq.tpr \
    -n ./build/index.ndx \
    -sl "$NZ" \
    -b "$BE_TIME" \
    -o ./build/density.xvg >> ./result/b_model.log 2>&1; then
    echo -e "${OK}Density calculation completed" | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}gmx density execution failed${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

# =============================
# Step 4: Extract average density in Z range
# =============================
echo "[STEP 6]Analyzing density in Z range [$BE_Z - $END_Z] nm" | tee -a ./result/b_model.log >&2

density=$(gmx analyze -f ./build/density.xvg -b "$BE_Z" -e "$END_Z" 2>/dev/null | grep "^SS1" | awk '{print $2}')

# Check result
if [[ "$density" == "NaN" ]]; then
    echo -e "${YELLOW}[WARNING]No valid data found in Z range [$BE_Z, $END_Z]${NC}" | tee -a ./result/b_model.log >&2
    density="undefined"
else
    echo -e "${OK}[RESULT]Average density = $density  kg/m³" | tee -a ./result/b_model.log >&2
    echo "OUTPUT: $density"
fi

# Save result to file
echo "$density" >> ./result/density_result.dat
echo -e "${GREEN}[FINISH]Pre-equilibration and density analysis completed${NC}" | tee -a ./result/b_model.log >&2