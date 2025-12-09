#!/bin/bash

source "${build_gmx}"

# Script Function: Calculate dipole angles
# Usage: script_name start_num end_num
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)


num=${analyze_begin_case}
end_num=${analyze_end_case}
count=0
execu_bin="${bash_dir}/../../bin/${analyze_cpp}"
# execu_bin="gmx density"
# export analyze_cpp="gmx_density"
src_file="${bash_dir}/module_C++/${analyze_cpp}.cpp"
mkdir -p ./deal_data/${analyze_cpp}
bash ${bash_dir}/../build_gmx_cpp.sh "${src_file}" "${execu_bin}"

# Process each case in loop



while (( "$num" <= "$end_num" )); do
    echo "Processing case$num..."
    if [[ ${isRerun} -eq 1 ]]; then
        PP_command="-f "./case$num/rerun_case/${DEFFNM}.${xtcORtrr}" -s "./case$num/rerun_case/${DEFFNM}.tpr" -n "./case$num/rerun_case/index.ndx""
    else
        PP_command="-f "./case$num/${DEFFNM}.${xtcORtrr}" -s "./case$num/${DEFFNM}.tpr" -n "./case$num/index.ndx""
        # PP_command="-f "./case$num/frame${num}.${xtcORtrr}" -s "./case$num/${DEFFNM}.tpr" -n "./case$num/index.ndx""
    fi
    
    out_command="-o "./deal_data/${analyze_cpp}/${num}${analyze_cpp}.xvg""
    if echo -e ${analyze_mol} | ${execu_bin}  ${PP_command} ${out_command} ${analysis_extra_command} -e 1 >> ./result/analyze.log 2>&1; then
        echo -e "${OK} command works fine. Starting analysis..."
    else
        echo -e "${YELLOW}WARNING! Incorrect location—this may be due to an error in the molecule name or duplicate molecule names in the index.${NC}"
        echo -e "Try fixing the duplicate molecule names."
        bash ${bash_dir}/../fix_duplicate_molname.sh ${num}
        echo -e "Retrying analysis..."
        if echo  ${analyze_mol} | ${execu_bin}  ${PP_command} ${out_command} ${analysis_extra_command} -e 1 >> ./result/analyze.log 2>&1; then
            echo -e "${OK} command works fine. Starting analysis..."
        else 
            echo -e "${ERROR} command failed in case${num} again!${NC}"
            # exit 1
            ((num++))
            continue;
        fi
    fi

    echo ${execu_bin}  ${PP_command} ${out_command} ${analysis_extra_command} >> ./result/analyze.log 2>&1
    {
        echo  ${analyze_mol} | ${execu_bin}  ${PP_command} ${out_command} ${analysis_extra_command}
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
rm -rf ./deal_data/${analyze_cpp}/*#
g++ "${bash_dir}/module_C++/average_xvg.cpp" -o "${bash_dir}/../../bin/ave_xvg" -O3 -std=c++17

if ${bash_dir}/../../bin/ave_xvg | tee debug ; then
    echo -e "${GREEN}Averaging xvg completed successfully!${NC}"
else
    echo "error in averaging xvg"
    exit 1
fi

for (( i=${analyze_begin_case}; i<=${analyze_end_case}; i++ )); do
rm ./deal_data/${analyze_cpp}/${i}${analyze_cpp}.xvg
done

echo -e "${GREEN}All tasks completed!${NC}"