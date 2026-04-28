#!/bin/bash

set -euo pipefail

other_condition() {
    if (( onlyUP == 1 )); then
        [ $(echo "$density <= $set_density" | bc -l) -eq 1 ]
    else
        false
    fi
}

get_molecule_counts() {
    if [[ ! -f "${TOP}.top" ]]; then
        echo "TOPOLOGY_NOT_FOUND"
        return
    fi
    awk '
        BEGIN { in_molecules=0 }
        {
            line=$0
            sub(/;.*/, "", line)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
            if (line == "") next
            if (line ~ /^\[/) {
                low=line
                gsub(/[[:space:]]/, "", low)
                low=tolower(low)
                in_molecules=(low=="[molecules]")
                next
            }
            if (in_molecules) {
                n=split(line, arr, /[[:space:]]+/)
                if (n >= 2) {
                    printf "%s=%s ", arr[1], arr[2]
                }
            }
        }
        END { print "" }
    ' "${TOP}.top"
}

log_build_model_info() {
    local iter_id="$1"
    local density_value="$2"
    local molecule_counts
    molecule_counts=$(get_molecule_counts)
    if [[ -z "$molecule_counts" ]]; then
        molecule_counts="NO_MOLECULE_INFO"
    fi
    printf "iter=%s density=%s %s\n" "$iter_id" "$density_value" "$molecule_counts" >> ./result/build_model_info.dat
}

SCRIPT_DIR="$(builtin cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
abs_error=$(echo "scale=5; $set_density * $error / 100" | bc -l)
build_dir=${SCRIPT_DIR}/../bin/build_cpp.sh
remake_src_cpp=${SCRIPT_DIR}/source/remake_packmol.cpp
remake_execu_bin=${SCRIPT_DIR}/../bin/remake_packmol
delete_src_cpp=${SCRIPT_DIR}/source/delete_molecules.cpp
delete_execu_bin=${SCRIPT_DIR}/../bin/delete_molecules
{
echo -e "${BLUE}"
echo "**********************************************************************"
echo "     I.   build model"
echo "⚙️  Parameters: target = $set_density ± $abs_error, max_iter = $max_iter"
echo "**********************************************************************"
echo -e "${NC}"
} | tee -a ./result/b_model.log >&2
# intial iter
iter=0
export iter
backup_top_file=".top"
backup_packmol_file=".inp"

if [[ ! -f "$backup_top_file" ]]; then
    cp "${TOP}.top" "$backup_top_file"
fi
if [[ ! -f "$backup_packmol_file" ]]; then
    cp "${packmol}.inp" "$backup_packmol_file"
fi
{
    echo "# build_model iteration info"
    echo "# format: iter=<n> density=<value> mol1=<count> mol2=<count> ..."
} > ./result/build_model_info.dat
# begin model
echo "Initializing system with packmol..." | tee -a ./result/b_model.log >&2
name_value=$(bash "$SCRIPT_DIR"/source/begin_packmol.sh | sed -n 's/^Final structure: \(.*\)\.gro$/\1/p' )
if [ -z "$name_value" ]; then
    echo -e "${ERROR}Error: Failed to run begin_packmol.sh or no valid output${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi
export name_value

# get initial density
output=$(sh "$SCRIPT_DIR"/source/check_density.sh)
if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
    echo -e "${ERROR}Error: Failed to get density from check_density.sh${NC}" | tee -a ./result/b_model.log >&2
    echo "Output was: $output" | tee -a ./result/b_model.log
    exit 1
fi

density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
    echo -e "${ERROR}Error: Invalid density value returned: $density${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi
density=$(printf "%.10f" "$density")
echo "Initial density: $density" | tee -a ./result/b_model.log >&2
log_build_model_info "$iter" "$density"
# main iter circle
while [ $iter -lt $max_iter ] && { (( $(echo "define abs(x) { if (x < 0) return -x; return x; }; abs($density - $set_density) > $abs_error" | bc -l) )) || other_condition; }; do
    iter=$((iter + 1))
    echo "Iteration $iter: Current density = $density, Target = $set_density" | tee -a ./result/b_model.log >&2
    # set factor
    diff=$(echo "$density - $set_density" | bc -l)
    scale_factor=$(echo "$set_density / $density" | bc -l)
    export scale_factor

    if (( $(echo "$diff > 0" | bc -l) )); then
        echo "Density too high, removing molecules directly with factor $scale_factor"  | tee -a ./result/b_model.log >&2
        if [ ! -f "${delete_execu_bin}" ]; then
            echo "Executable ${delete_execu_bin} does not exist, starting compilation..."
            bash "${build_dir}" "${delete_src_cpp}" "${delete_execu_bin}"
            if [ $? -eq 0 ] && [ -f "${delete_execu_bin}" ]; then
                echo "Compilation successful!" | tee -a ./result/b_model.log >&2
            else
                echo -e "${ERROR}Compilation failed!" | tee -a ./result/b_model.log >&2
                exit 1
            fi
        else
            echo "Executable ${delete_execu_bin} already exists, skipping compilation." | tee -a ./result/b_model.log >&2
        fi

        if ! "${delete_execu_bin}" "$scale_factor" "${name_value}.gro" "${TOP}.top" "${packmol}.inp" | tee -a ./result/b_model.log >&2; then
            echo -e "${ERROR}Error: delete_molecules failed with scale factor $scale_factor${NC}" | tee -a ./result/b_model.log >&2
            exit 1
        fi
    else
        echo "Density too low, rebuilding system with factor $scale_factor" | tee -a ./result/b_model.log >&2
        if [ ! -f "${remake_execu_bin}" ]; then
            echo "Executable ${remake_execu_bin} does not exist, starting compilation..."
            bash "${build_dir}" "${remake_src_cpp}" "${remake_execu_bin}"
            if [ $? -eq 0 ] && [ -f "${remake_execu_bin}" ]; then
                echo "Compilation successful!" | tee -a ./result/b_model.log >&2
            else
                echo -e "${ERROR}Compilation failed!" | tee -a ./result/b_model.log >&2
                exit 1
            fi
        else
            echo "Executable ${remake_execu_bin} already exists, skipping compilation." | tee -a ./result/b_model.log >&2
        fi

        if  ! "${remake_execu_bin}" "$scale_factor" "${packmol}" | tee -a ./result/b_model.log >&2 ; then
            echo -e "${ERROR}Error: remake_packmol failed with scale factor $scale_factor${NC}"  | tee -a ./result/b_model.log >&2
            exit 1
        fi

        new_name_value=$(sh "$SCRIPT_DIR"/source/begin_packmol.sh | sed -n 's/^Final structure: \(.*\)\.gro$/\1/p')
        if [ -z "$new_name_value" ]; then
            echo -e "${ERROR}Error: begin_packmol.sh finished but no valid structure name was returned${NC}" | tee -a ./result/b_model.log >&2
            exit 1
        fi
        name_value="$new_name_value"
        export name_value
    fi

    # calculate density
    output=$(sh "$SCRIPT_DIR"/source/check_density.sh)
    if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
        echo -e "${ERROR}Error: Failed to get updated density.${NC}" | tee -a ./result/b_model.log  >&2
        echo "Output was: $output"  | tee -a ./result/b_model.log 
        exit 1
    fi
    density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
    if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
        echo -e "${ERROR}Error: Invalid density value returned: $density${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi
    density=$(printf "%.10f" "$density")
    log_build_model_info "$iter" "$density"
done

# result
if [ $iter -ge $max_iter ]; then
    echo -e "${ERROR}error: Maximum iterations ($max_iter) reached. Density did not converge.${NC}"  | tee -a ./result/b_model.log >&2
    echo "Final density: $density; but target is $set_density ± $abs_error" | tee -a ./result/b_model.log  >&2
    exit 1
else
    echo -e "${GREEN}Success: Density converged to $density within $iter iterations.${NC}" | tee -a ./result/b_model.log  >&2
    echo -e "Target: ${set_density} ± ${abs_error}" >&2
fi
rm "${name_value}".gro
