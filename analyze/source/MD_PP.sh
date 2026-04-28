#!/bin/bash

set -euo pipefail

if [[ -n "${build_gmx:-}" && -f "${build_gmx}" ]]; then
    source "${build_gmx}"
fi

# Script Function: Calculate dipole angles
# Usage: script_name start_num end_num
bash_dir=$(builtin cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)


num=${analyze_begin_case}
end_num=${analyze_end_case}
count=0
analysis_program_name="${analyze_cpp}"
analysis_program_label="${analysis_program_name//[^[:alnum:]_.-]/_}"
src_file="${bash_dir}/module_C++/${analysis_program_name}.cpp"

declare -a execu_cmd extra_cmd
if [[ -f "${src_file}" ]]; then
    execu_bin="${bash_dir}/../../bin/${analysis_program_name}"
    execu_cmd=("${execu_bin}")
    bash "${bash_dir}/../build_gmx_cpp.sh" "${src_file}" "${execu_bin}"
else
    read -r -a execu_cmd <<< "${analysis_program_name}"
    if (( ${#execu_cmd[@]} == 0 )); then
        echo -e "${ERROR} Empty analysis command.${NC}"
        exit 1
    fi
    if ! command -v "${execu_cmd[0]}" >/dev/null 2>&1; then
        echo -e "${ERROR} External analysis command '${execu_cmd[0]}' not found.${NC}"
        exit 1
    fi
fi

read -r -a extra_cmd <<< "${analysis_extra_command:-}"
export analyze_cpp="${analysis_program_label}"
mkdir -p "./deal_data/${analysis_program_label}"

run_pp_command() {
    local case_num="$1"
    local traj_dir
    local output_file="./deal_data/${analysis_program_label}/${case_num}${analysis_program_label}.xvg"
    local -a cmd

    if [[ ${isRerun} -eq 1 ]]; then
        traj_dir="./case${case_num}/rerun_case"
    else
        traj_dir="./case${case_num}"
    fi

    cmd=(
        "${execu_cmd[@]}"
        -f "${traj_dir}/${DEFFNM_analyze}.${xtcORtrr}"
        -s "${traj_dir}/${DEFFNM_analyze}.tpr"
        -n "${traj_dir}/index.ndx"
        -o "${output_file}"
        "${extra_cmd[@]}"
        -e 1
    )

    printf '%b\n' "${analyze_mol}" | "${cmd[@]}"
}

log_pp_command() {
    local case_num="$1"
    local traj_dir
    local output_file="./deal_data/${analysis_program_label}/${case_num}${analysis_program_label}.xvg"
    local -a cmd

    if [[ ${isRerun} -eq 1 ]]; then
        traj_dir="./case${case_num}/rerun_case"
    else
        traj_dir="./case${case_num}"
    fi

    cmd=(
        "${execu_cmd[@]}"
        -f "${traj_dir}/${DEFFNM_analyze}.${xtcORtrr}"
        -s "${traj_dir}/${DEFFNM_analyze}.tpr"
        -n "${traj_dir}/index.ndx"
        -o "${output_file}"
        "${extra_cmd[@]}"
        -e 1
    )

    printf '%q ' "${cmd[@]}" >> ./result/analyze.log
    printf '\n' >> ./result/analyze.log
}

# Process each case in loop



while (( "$num" <= "$end_num" )); do
    echo "Processing case$num..."
    if run_pp_command "${num}" >> ./result/analyze.log 2>&1; then
        echo -e "${OK} command works fine. Starting analysis..."
    else
        echo -e "${YELLOW}WARNING! Incorrect location—this may be due to an error in the molecule name or duplicate molecule names in the index.${NC}"
        echo -e "Try fixing the duplicate molecule names."
        bash "${bash_dir}/../fix_duplicate_molname.sh" "${num}"
        echo -e "Retrying analysis..."
        if run_pp_command "${num}" >> ./result/analyze.log 2>&1; then
            echo -e "${OK} command works fine. Starting analysis..."
        else 
            echo -e "${ERROR} command failed in case${num} again!${NC}"
            # exit 1
            ((num++))
            continue;
        fi
    fi

    log_pp_command "${num}"
    {
        run_pp_command "${num}"
    } >> ./result/analyze.log 2>&1 &
    #     echo "Processing case$num..."
    # else
    #     echo -e "${ERROR} command failed in case${num}!${NC}"
    #     exit 1
    # fi
    # Simple task control: wait after every analyze_core task
    ((count++))
    if (( "$count" == "${analyze_core}" )) ; then
        wait
        count=0
        echo -e "${GREEN}Processed up to case$num. Continuing...${NC}"
    fi
    ((num++))
done

# Wait for all background tasks to complete
wait
rm -rf "./deal_data/${analysis_program_label}"/*#
g++ "${bash_dir}/module_C++/average_xvg.cpp" -o "${bash_dir}/../../bin/ave_xvg" -O3 -std=c++17

if ${bash_dir}/../../bin/ave_xvg | tee debug ; then
    echo -e "${GREEN}Averaging xvg completed successfully!${NC}"
else
    echo "error in averaging xvg"
    exit 1
fi

for (( i=${analyze_begin_case}; i<=${analyze_end_case}; i++ )); do
rm "./deal_data/${analysis_program_label}/${i}${analysis_program_label}.xvg"
done

echo -e "${GREEN}All tasks completed!${NC}"
