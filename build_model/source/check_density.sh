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
# Step 2: mdrun - run
# =============================
echo "[STEP 2]Running mdrun for pre-equilibration" | tee -a ./result/b_model.log >&2
mini_pre_eq_overlap_risk=0

set +o pipefail
stdbuf -o0 gmx mdrun \
    -s ./build/mini_pre_eq.tpr \
    -deffnm ./build/mini_pre_eq \
    -ntmpi 1 \
    -ntomp "$NPOS" \
    -v 2>&1 | stdbuf -o0 tee -a ./result/b_model.log | awk 'BEGIN{RS="\r|\n"} /^step/{printf "\r\033[K%s", $0 > "/dev/stderr"; fflush("/dev/stderr")} /^Performance/{printf "\n%s\n", $0 > "/dev/stderr"; fflush("/dev/stderr")}'
mdrun_rc=${PIPESTATUS[0]}
set -o pipefail
if [[ $mdrun_rc -ne 0 ]]; then
    echo -e "${ERROR}min energy failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi
echo -e "${OK}min energy completed successfully" | tee -a ./result/b_model.log >&2
mini_pre_eq_log=./build/mini_pre_eq.log
if [[ -f "$mini_pre_eq_log" ]]; then
    if grep -q "did not reach the requested Fmax" "$mini_pre_eq_log"; then
        mini_pre_eq_overlap_risk=1
        echo -e "${YELLOW}[WARNING]mini_pre_eq reached machine precision but did not reach requested Fmax.${NC}" | tee -a ./result/b_model.log >&2
        echo -e "${YELLOW}[WARNING]Continue workflow. If next MD fails with 'The total potential energy is ... extremely high', script will stop and report overlap/box issue.${NC}" | tee -a ./result/b_model.log >&2
    fi
else
    echo -e "${YELLOW}[WARNING]mini_pre_eq log not found: ${mini_pre_eq_log}. Skip Fmax safeguard check.${NC}" | tee -a ./result/b_model.log >&2
fi

# Rename output file
mv ./build/mini_pre_eq.gro ./build/pre_eq.gro


# =============================
# Step 1: grompp - generate tpr
# =============================
echo "[STEP 3]Running grompp to generate pre_eq.tpr" | tee -a ./result/b_model.log >&2

cp "${MDP}".mdp ./build/grompp.mdp
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_den}/" ./build/grompp.mdp

eq_input_gro="./build/pre_eq.gro"
if [[ "${ENABLE_ANNEAL:-no}" == "yes" ]]; then
    echo "[STEP 3A]Annealing enabled. Preparing anneal.mdp (T=${ANNEAL_TEMP:-400} K)" | tee -a ./result/b_model.log >&2
    cp ./build/grompp.mdp ./build/anneal.mdp

    if [[ -n "${ANNEAL_NSTEPS:-}" ]]; then
        anneal_nsteps="${ANNEAL_NSTEPS}"
    else
        echo -e "${ERROR}ANNEAL_NSTEPS not set. "| tee -a ./result/b_model.log >&2
        exit 1
    fi

    if ! [[ "$anneal_nsteps" =~ ^[0-9]+$ ]] || [[ "$anneal_nsteps" -le 0 ]]; then
        echo -e "${ERROR}Invalid anneal nsteps: ${anneal_nsteps}. Check ANNEAL_TIME_PS/ANNEAL_NSTEPS.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    if grep -Eq "^[[:space:]]*ref[-_]t[[:space:]]*=" ./build/anneal.mdp; then
        sed -i -E "s/^[[:space:]]*ref[-_]t[[:space:]]*=.*/ref_t = ${ANNEAL_TEMP:-400}/" ./build/anneal.mdp
    else
        echo -e "${ERROR}invalid anneal.mdp: ref_t parameter not found. " | tee -a ./result/b_model.log >&2
        exit 1
    fi

    sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${anneal_nsteps}/" ./build/anneal.mdp

    echo "[STEP 3B]Running grompp for annealing (nsteps=${anneal_nsteps})" | tee -a ./result/b_model.log >&2
    if gmx grompp \
        -f ./build/anneal.mdp \
        -c "$eq_input_gro" \
        -p "${TOP}.top" \
        -n ./build/index.ndx \
        -o ./build/anneal.tpr \
        -maxwarn "${maxWarn}" >> ./result/b_model.log 2>&1; then
        echo -e "${OK}anneal grompp completed successfully" | tee -a ./result/b_model.log >&2
    else
        echo -e "${ERROR}anneal grompp failed.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    echo "[STEP 3C]Running mdrun for annealing" | tee -a ./result/b_model.log >&2
    if [[ "$GPU" -eq 1 ]]; then
        echo -e "\tUsing GPU acceleration for annealing" | tee -a ./result/b_model.log >&2
        anneal_mdrun_cmd=(
            gmx mdrun
            -s ./build/anneal.tpr
            -deffnm ./build/anneal
            -ntmpi 1
            -ntomp "$NPOS"
            -pme gpu
            -pmefft gpu
            -nb gpu
        )
        [[ "${ENABLE_BONDED_GPU:-1}" -eq 1 ]] && anneal_mdrun_cmd+=(-bonded gpu)
        anneal_mdrun_cmd+=(
            -tunepme no
            -v
        )
    else
        echo -e "\tUsing CPU computation for annealing" | tee -a ./result/b_model.log >&2
        anneal_mdrun_cmd=(
            gmx mdrun
            -s ./build/anneal.tpr
            -deffnm ./build/anneal
            -ntmpi 1
            -ntomp "$NPOS"
            -v
        )
    fi
    
    set +o pipefail
    stdbuf -o0 "${anneal_mdrun_cmd[@]}" 2>&1 | stdbuf -o0 tee -a ./result/b_model.log | awk 'BEGIN{RS="\r|\n"} /^step/{printf "\r\033[K%s", $0 > "/dev/stderr"; fflush("/dev/stderr")} /^Performance/{printf "\n%s\n", $0 > "/dev/stderr"; fflush("/dev/stderr")}'
    mdrun_rc=${PIPESTATUS[0]}
    set -o pipefail
    if [[ $mdrun_rc -ne 0 ]]; then
        echo -e "${ERROR}anneal mdrun failed (exit code: $mdrun_rc)${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi

    if [[ -s ./build/anneal.gro ]]; then
        echo -e "${OK}annealing completed successfully, anneal.gro will be used for normal pre-equilibration" | tee -a ./result/b_model.log >&2
        eq_input_gro="./build/anneal.gro"
    else
        echo -e "${ERROR}annealing failed: anneal.gro not generated or empty${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi
else
    echo "[STEP 3A]Annealing disabled (ENABLE_ANNEAL=${ENABLE_ANNEAL:-no}). Continue with normal pre-equilibration." | tee -a ./result/b_model.log >&2
fi

if gmx grompp \
    -f ./build/grompp.mdp \
    -c "$eq_input_gro" \
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
    )
    [[ "${ENABLE_BONDED_GPU:-1}" -eq 1 ]] && mdrun_cmd+=(-bonded gpu)
    mdrun_cmd+=(
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
mdrun_failed=0

set +o pipefail
stdbuf -o0 "${mdrun_cmd[@]}" 2>&1 | stdbuf -o0 tee -a ./result/b_model.log | awk 'BEGIN{RS="\r|\n"} /^step/{printf "\r\033[K%s", $0 > "/dev/stderr"; fflush("/dev/stderr")} /^Performance/{printf "\n%s\n", $0 > "/dev/stderr"; fflush("/dev/stderr")}'
mdrun_rc=${PIPESTATUS[0]}
set -o pipefail
if [[ $mdrun_rc -ne 0 ]]; then
    echo -e "${YELLOW}[WARNING]mdrun exited with code: $mdrun_rc${NC}" | tee -a ./result/b_model.log >&2
    mdrun_failed=1
fi

# Execute mdrun
if [[ "$mdrun_failed" -eq 0 ]] && [[ -s ./build/pre_eq.xtc ]]; then
    echo -e "${OK}mdrun completed successfully" | tee -a ./result/b_model.log >&2
else
    pre_eq_log=./build/pre_eq.log
    if [[ -f "$pre_eq_log" ]] \
        && grep -q "Fatal error:" "$pre_eq_log" \
        && grep -Eq "The total potential energy is .*extremely high" "$pre_eq_log"; then
        if [[ "$mini_pre_eq_overlap_risk" -eq 1 ]]; then
            echo -e "${YELLOW}[WARNING]mini_pre_eq.log already reported 'did not reach the requested Fmax', indicating strong initial overlap risk.${NC}" | tee -a ./result/b_model.log >&2
            echo -e "${ERROR} Energy minimization failed. Your molecules are too much OR Your box is too small!${NC}" | tee -a ./result/b_model.log >&2
        fi
    else
        echo -e "${ERROR}mdrun failed, Unknown error${NC}" | tee -a ./result/b_model.log >&2
    fi
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
