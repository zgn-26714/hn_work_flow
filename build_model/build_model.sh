#!/bin/bash


other_condition() {
    if (( onlyUP == 1 )); then
        [ $(echo "$density <= $set_density" | bc -l) -eq 1 ]
    else
        true
    fi
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
abs_error=$(echo "scale=5; $set_density * $error / 100" | bc -l)
{
echo -e "${BLUE}"
echo "**********************************************************************"
echo "     I.   build model"
echo "⚙️  Parameters: target = $set_density ± $abs_error, max_iter = $max_iter, factor = $factor"
echo "**********************************************************************"
echo -e "${NC}"
} | tee -a ./result/b_model.log >&2
# intial iter
iter=0
export iter
cp "${packmol}".inp packmol_.inp
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

cp "${TOP}".top baktop1.top
# main iter circle
while { (( $(echo "define abs(x) { if (x < 0) return -x; return x; }; abs($density - $set_density) > $abs_error" | bc -l) )) && [ $iter -lt $max_iter ]; } && other_condition; do
    iter=$((iter + 1))
    echo "Iteration $iter: Current density = $density, Target = $set_density" | tee -a ./result/b_model.log >&2
    # set factor
    diff=$(echo "$density - $set_density" | bc -l)
    if (( $(echo "$diff > 0" | bc -l) )); then
        scale_factor=$(echo "$set_density / $density" | bc -l)
        echo "Density too high, expanding system by factor $scale_factor"  | tee -a ./result/b_model.log >&2
    else
        scale_factor=$(echo "$density / $set_density" | bc -l)
        echo "Density too low, compressing system by factor $scale_factor" | tee -a ./result/b_model.log >&2
    fi
    #remake packmol and top
    export scale_factor
    build_dir=${SCRIPT_DIR}/../bin/build_cpp.sh
    src_cpp=${SCRIPT_DIR}/source/remake_packmol.cpp
    execu_bin=${SCRIPT_DIR}/../bin/remake_packmol
    if [ ! -f "${execu_bin}" ]; then
        echo "Executable ${execu_bin} does not exist, starting compilation..."
        bash ${build_dir} ${src_cpp} ${execu_bin}
        
        # Check if compilation was successful
        if [ $? -eq 0 ] && [ -f "${execu_bin}" ]; then
            echo "Compilation successful!"
        else
            echo "Compilation failed!"
            exit 1
        fi
    else
        echo "Executable ${execu_bin} already exists, skipping compilation."
    fi
    if  ${execu_bin} "$scale_factor" packmol_; then
        echo -e "${ERROR}Error: remake_packmol failed with scale factor $scale_factor${NC}"  | tee -a ./result/b_model.log >&2
        exit 1
    fi
    sh "$SCRIPT_DIR"/source/begin_packmol.sh  | tee -a ./result/b_model.log >&2
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
mv packmol_.inp result.inp