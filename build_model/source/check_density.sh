#!/bin/bash

# =============================
# Safety settings and checks
# =============================
source "$GMXRC"
set -euo pipefail
echo "
===================================================================================================
 ---------------------------------------- iter="${iter}" -----------------------------------------
  "【INFO】Starting pre-equilibration and density analysis workflow"
  Script    : ${BASH_SOURCE[0]}
===================================================================================================">&2

# Check if GMXRC exists
if [[ ! -f "$GMXRC" ]]; then
    echo "【ERROR】GMXRC file not found: $GMXRC" >&2
    exit 1
fi


# Clean up old files safely
echo "【INFO】Cleaning up temporary files...">&2
rm -f pre_eq*
rm -f \#*




echo "【STEP 1】Running grompp to generate mini_pre_eq.tpr">&2

if gmx grompp \
    -f "${mini_MDP}.mdp" \
    -c "${name_value}.gro" \
    -p "${TOP}.top" \
    -n "${INDEX}.ndx" \
    -o mini_pre_eq.tpr \
    -maxwarn "${maxWarn}" > test.log 2>&1; then
    echo "【SUCCESS】grompp completed successfully">&2
else
    echo "【ERROR】grompp failed.">&2
    exit 1
fi

# =============================
# Step 2: mdrun - run simulation
# =============================
echo "【STEP 2】Running mdrun for pre-equilibration">&2
gmx mdrun \
    -s mini_pre_eq.tpr \
    -deffnm mini_pre_eq \
    -ntmpi 1 \
    -ntomp "$NPOS" \
    -v > test.log 2>&1

# Rename output file
mv mini_pre_eq.gro pre_eq.gro


# =============================
# Step 1: grompp - generate tpr
# =============================
echo "【STEP 3】Running grompp to generate pre_eq.tpr">&2
if gmx grompp \
    -f "${MDP}.mdp" \
    -c pre_eq.gro \
    -p "${TOP}.top" \
    -n "${INDEX}.ndx" \
    -o pre_eq.tpr \
    -maxwarn "${maxWarn}" > test.log 2>&1; then
    echo "【SUCCESS】grompp completed successfully">&2
else
    echo "【ERROR】grompp failed.">&2
    exit 1
fi

# =============================
# Step 2: mdrun - run simulation
# =============================
echo "【STEP 4】Running mdrun for pre-equilibration">&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\t【INFO】Using GPU acceleration">&2
    mdrun_cmd=(
        gmx mdrun
        -s pre_eq.tpr
        -deffnm pre_eq
        -ntmpi 1
        -ntomp "$NPOS"
        -pme gpu
        -pmefft gpu
        -nb gpu
        -tunepme no
        -v
    )
else
    echo -e "\t【INFO】Using CPU computation">&2
    mdrun_cmd=(
        gmx mdrun
        -s pre_eq.tpr
        -deffnm pre_eq
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi
"${mdrun_cmd[@]}" > test.log 2>&1 

# Execute mdrun
if [[ -s pre_eq.xtc ]]; then
    echo "【SUCCESS】mdrun completed successfully">&2
else
    echo "【ERROR】mdrun failed: pre_eq.edr not generated or empty">&2
    exit 1
fi

# =============================
# Step 3: Compute density profile
# =============================
echo "【STEP 5】Computing density profile">&2
if echo 0|gmx density \
    -f pre_eq.xtc \
    -s pre_eq.tpr \
    -n "${INDEX}.ndx" \
    -sl "$NZ" \
    -b "$BE_TIME" \
    -o density.xvg > test.log 2>&1; then
    echo "【SUCCESS】Density calculation completed">&2
else
    echo "【ERROR】gmx density execution failed">&2
    exit 1
fi

# =============================
# Step 4: Extract average density in Z range
# =============================
echo "【STEP 6】Analyzing density in Z range [$BE_Z - $END_Z] nm">&2

density=$(gmx analyze -f density.xvg -b "$BE_Z" -e "$END_Z" 2>/dev/null | grep "^SS1" | awk '{print $2}')

# Check result
if [[ "$density" == "NaN" ]]; then
    echo "【WARNING】No valid data found in Z range [$BE_Z, $END_Z]">&2
    density="undefined"
else
    echo "【RESULT】Average density = $density  kg/m³ ">&2
    echo "OUTPUT: $density"
fi

# Save result to file
echo "$density" >> density_result.txt
echo "【FINISH】Pre-equilibration and density analysis completed">&2