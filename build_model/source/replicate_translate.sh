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