#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
abs_error=$(echo "$set_density * $error / 100" | bc -l)
echo ""
echo "**********************************************************************"
echo "     I.   build model"
echo "⚙️  Parameters: target = $set_density ± $abs_error, max_iter = $max_iter, factor = $factor"
echo "**********************************************************************"
echo ""
# intial iter
iter=0
export iter
# begin model
echo "Initializing system with packmol..."
name_value=$(sh "$SCRIPT_DIR"/source/begin_packmol.sh | sed -n 's/^Final structure: \(.*\)\.gro$/\1/p' )
if [ -z "$name_value" ]; then
    echo "Error: Failed to run begin_packmol.sh or no valid output" >&2
    exit 1
fi
export name_value
rm -f density_result.dat

# get initial density
output=$(sh "$SCRIPT_DIR"/source/check_density.sh)
if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
    echo "Error: Failed to get density from check_density.sh" >&2
    echo "Output was: $output"
    exit 1
fi

density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
    echo "Error: Invalid density value returned: $density" >&2
    exit 1
fi
density=$(printf "%.10f" "$density")
echo "Initial density: $density"

# main iter circle
while (( $(echo "define abs(x) { if (x < 0) return -x; return x; }; abs($density - $set_density) > $abs_error" | bc -l) )) && [ $iter -lt $max_iter ]; do
    iter=$((iter + 1))
    echo "Iteration $iter: Current density = $density, Target = $set_density" >&2
    # set factor
    diff=$(echo "$density - $set_density" | bc -l)
    if (( $(echo "$diff > 0" | bc -l) )); then
        scale_factor=$(echo "1 - $factor" | bc -l)
        echo "Density too high, expanding system by factor $scale_factor"
    else
        scale_factor=$(echo "1 + $factor" | bc -l)
        echo "Density too low, compressing system by factor $scale_factor"
    fi
    # remake_packmol.py
    export scale_factor
    if ! python3 "$SCRIPT_DIR"/source/remake_packmol.py "$scale_factor"; then
        echo "Error: remake_packmol.py failed with scale factor $scale_factor" >&2
        exit 1
    fi
    sh "$SCRIPT_DIR"/source/begin_packmol.sh >&2
    # calculate density
    output=$(sh "$SCRIPT_DIR"/source/check_density.sh)
    if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
        echo "Error: Failed to get updated density." >&2
        echo "Output was: $output"
        exit 1
    fi
    density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
    if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
        echo "Error: Invalid density value returned: $density" >&2
        exit 1
    fi
    density=$(printf "%.10f" "$density")
done

# result
if [ $iter -ge $max_iter ]; then
    echo "Warning: Maximum iterations ($max_iter) reached. Density did not converge." >&2
    echo "Final density: $density" >&2
    exit 1
else
    echo "Success: Density converged to $density within $iter iterations." >&2
    echo "Target: $set_density ± $error" >&2
fi