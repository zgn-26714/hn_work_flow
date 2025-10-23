#!/bin/bash

# Script Function: Calculate dipole angles
# Usage: script_name start_num end_num
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

num=${analyze_begin_case}
end_num=${analyze_end_case}
PP_command="-f "./case$num/${DEFFNM}.xtc" -s "./case$num/${DEFFNM}.tpr" -n "./case$num/index.ndx""
out_command="-o "./deal_data/dipole_ang/${num}angle_dipole_z.xvg""
count=0

execu_bin="${bash_dir}/../../bin/${analyze_cpp}"
src_file="${bash_dir}/module_C++/${analyze_cpp}.cpp"
mkdir -p ./deal_data/${analyze_cpp}
bash ${bash_dir}/../build_cpp.sh "${src_file}" "${execu_bin}"

# Process each case in loop
while (( "$num" <= "$end_num" )); do
    if echo  ${analyze_mol} | ${execu_bin}  ${PP_command} ${out_command} ${analysis_extra_command} & ; then;
        echo "Processing case$num..."
    else
        echo -e "${ERROR} command failed in case${num}!${NC}"
    fi
    # Simple task control: wait after every analyze_core task
    ((count++))
    if [ "$count" -eq ${analyze_core} ]; then
        wait
        count=0
    fi
    ((num++))
done

# Wait for all background tasks to complete
wait
echo -e "${GREEN}All tasks completed!${NC}"