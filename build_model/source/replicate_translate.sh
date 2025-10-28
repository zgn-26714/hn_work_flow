#!/bin/bash

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if gmx editconf -f ./build/pre_eq.gro -o ./build/pre_eq2.gro -translate 0 0 ${zbox} ;then
    echo -e "${OK}Successfully translated structure along z-axis by ${zbox} nm." | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}Translation failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

cpp_src=${bash_dir}/merge_gro.cpp
execu_bin=${bash_dir}/../../bin/replicate_translate
bash ${bash_dir}/build_cpp.sh ${cpp_src} ${execu_bin}

${execu_bin} "./build/pre_eq.gro" "./build/pre_eq2.gro"

formatted_line=$(printf "%8.4f%8.4f%8.4f" $xbox $ybox $zbox)
sed -i "\$ s/^.*$/      $formatted_line/" ./build/pre_merge.gro
echo -e "${OK}${GREEN}merge two gro successfully!${NC}" | tee -a ./result/b_model.log >&2

cpp_src=${bash_dir}/double_slit_top.cpp
execu_bin=${bash_dir}/../../bin/double_slit_top
bash ${bash_dir}/build_cpp.sh ${cpp_src} ${execu_bin}
${execu_bin} >> ./result/b_model.log 2>&1

#改名这里还需要很多的测试；主要是手动也没测试几次
if [ "$is_change_name" -eq 1 ]; then
    if [ -n "$need_change_names" ] && [ -n "$target_names" ]; then
        cpp_src=${bash_dir}/change_name.cpp
        execu_bin=${bash_dir}/../../bin/change_name
        bash ${bash_dir}/build_cpp.sh ${cpp_src} ${execu_bin}
        ${execu_bin} "./build/pre_merge.gro"  >> ./result/b_model.log 2>&1
    else
        echo -e "${ERROR}need_change_names or target_names not defined.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi
fi


