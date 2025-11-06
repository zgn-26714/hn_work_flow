#!/bin/bash

set -euo pipefail  # Strict error handling

# Simple calculator function
calc() {
    python -c "print($*)"
}

# ============== Configuration Parameters ===============
workdir="$(pwd)"
framedir="$workdir/${initialframes}"
matrixdir="$workdir/${matrixdata}"
basicdir="$workdir/${basicfile}"

# Validate required directories and files
validate_dependencies() {
    [[ -d "$framedir" ]] || { echo "❌ Frame directory not found: $framedir"; exit 1; }
    [[ -d "$matrixdir" ]] || { echo "❌ Matrix directory not found: $matrixdir"; exit 1; }
    [[ -d "$basicdir" ]] || { echo "❌ Basic file directory not found: $basicdir"; exit 1; }
}

# Set parameters with default values
num=${num}              # Number of cases to process per batch
queue=${queue}           # PBS queue name
cas=${START}              # Starting case number
endcas=${END}             # Ending case number

# Required parameters validation
if [[ -z "${Temperature:-}" || -z "${V:-}" || -z "${ic:-}" ]]; then
    echo "❌ Missing required parameters: Temperature, V, or ic" | tee -a ./result/run_md.log >&2
    exit 1
fi

# =============== Function Definitions ===============
make_case() {
    local case_num=$1
    local objdir="$workdir/charging/${Temperature}k/${V}V/${ic}ps/case${case_num}"
    local slowdir="${workdir}/slowdir"

    # Create directory and copy files
    mkdir -p "$objdir"
    cp "$framedir/frame${case_num}.gro" "$objdir/"
    cp "$basicdir"/* "$objdir/" 2>/dev/null || echo "⚠️  No files to copy from basicdir" | tee -a  ./result/run_md.log >&2
    cp "$matrixdir/CPM_ControlFile.dat_0V" "$objdir/CPM_ControlFile.dat"

    # Remove old files and create symbolic links
    rm -f "$objdir/allMatrixA.bin" "$objdir/Dphis_control.dat"
    ln -sf "$matrixdir/allMatrixA.bin" "$objdir/allMatrixA.bin"
    
    local dphis_file="$slowdir/Dphis_control.dat${V}_${ic}_${SKIPTIME_PS}"
    if [[ -f "$dphis_file" ]]; then
        ln -sf "$dphis_file" "$objdir/Dphis_control.dat"
    else
        echo "⚠️  Dphis file not found: $dphis_file" | tee -a  ./result/run_md.log >&2
    fi

    echo "Created case: $case_num at $objdir" | tee -a ./result/run_md.log >&2
}

# =============== Main Execution ===============
main() {
    echo "🚀 Starting case generation process..." | tee -a ./result/run_md.log >&2
    echo "📊 Parameters: Temperature=${Temperature}K, Voltage=${V}V, Time=${ic}ps" | tee -a ./result/run_md.log >&2
    echo "📁 Range: cases $cas to $endcas" | tee -a ./result/run_md.log >&2
    
    # Validate dependencies first
    validate_dependencies
    
    # Get script directory and run generator
    local bash_dir
    bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    
    if  bash "$bash_dir/source/generate_scanfile.sh" | tee -a ./result/run_md.log; then
        echo -e "${GREEN}✅Success generating scan file...${NC}" | tee -a ./result/run_md.log >&2
    else
        echo -e "${ERROR} generate_scanfile.sh not found or error, exit!" | tee -a ./result/run_md.log >&2
        eixt 1
    fi

    local case_dirs=()
    
    # Process cases in batches
    while (( cas <= endcas )); do
        echo "🔄 Processing batch starting from case $cas..." | tee -a ./result/run_md.log >&2   
        make_case "$cas"
        (( cas ++ )) 
    done

    # Summary
    echo -e "${GREEN}✅ Case generation completed successfully!${NC}" | tee -a ./result/run_md.log >&2
}

# Run main function
main "$@"