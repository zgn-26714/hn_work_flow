#!/bin/bash

set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if (( isSlit == 1)); then
    last_line=$(tail -n 1 ./model/single_electrode.gro)
    read -r xbox ybox zbox <<< "$last_line"
fi
if gmx editconf -f ./build/pre_eq.gro -o ./build/pre_eq2.gro -translate 0 0 ${zbox} >> ./result/b_model.log 2>&1;then
    echo -e "${OK}Successfully translated structure along z-axis by ${zbox} nm." | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}Translation failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

cpp_src=${bash_dir}/merge_gro.cpp
execu_bin=${bash_dir}/../../bin/replicate_translate
bash ${bash_dir}/../../bin/build_cpp.sh ${cpp_src} ${execu_bin}

${execu_bin} ./build/pre_eq.gro ./build/pre_eq2.gro

double_z_box=$(echo "$zbox * 2" | bc -l)
formatted_line=$(printf "%8.4f%8.4f%8.4f" $xbox $ybox $double_z_box)
sed -i "\$ s/^.*$/      $formatted_line/" ./build/pre_eq_merge.gro
echo -e "${OK}${GREEN}merge two gro successfully!${NC}" | tee -a ./result/b_model.log >&2

cpp_src=${bash_dir}/double_slit_top.cpp
execu_bin=${bash_dir}/../../bin/double_slit_top
bash ${bash_dir}/../../bin/build_cpp.sh ${cpp_src} ${execu_bin}
${execu_bin} >> ./result/b_model.log 2>&1
mv ${TOP}.top.processed ${TOP}.top
# #从generate_slit开始规范名字，关闭此功能
# if [ "$is_change_name" -eq 1 ]; then
#     if [ -n "$need_change_names" ] && [ -n "$target_names" ]; then
#         cpp_src=${bash_dir}/change_name.cpp
#         execu_bin=${bash_dir}/../../bin/change_name
#         bash ${bash_dir}/build_cpp.sh ${cpp_src} ${execu_bin}
#         ${execu_bin} "./build/pre_merge.gro"  >> ./result/b_model.log 2>&1
#     else
#         echo -e "${ERROR}need_change_names or target_names not defined.${NC}" | tee -a ./result/b_model.log >&2
#         exit 1
#     fi
# fi

echo "Rebalance..." | tee -a ./result/b_model.log >&2
rm -f *#
mv build/pre_eq_merge.gro build/pre_eq.gro

echo "[step 1] genenrate .tpr (mini_energy) (grompp)" | tee -a ./result/b_model.log >&2

cp "${mini_MDP}".mdp ./build/mini.mdp
sed -i -e 's/^freezegrps.*$/freezegrps = CL  CR  GRA/' -e 's/^freezedim.*$/freezedim   = Y Y Y Y Y Y Y Y Y/' ./build/mini.mdp

echo q| gmx make_ndx -f ./build/pre_eq.gro -o ./build/index.ndx >> ./result/b_model.log  2>&1


gmx grompp \
    -f ./build/mini.mdp \
    -c ./build/pre_eq.gro \
    -o ./build/mini_re.tpr \
    -maxwarn "${maxWarn}" \
    -p "${TOP}.top" \
    -n ./build/index.ndx  >> ./result/b_model.log  2>&1 \
    || { echo -e "${ERROR} gmx grompp failed (mini_energy)${NC}" | tee -a ./result/b_model.log >&2; exit 1; }

#  mdrun 
echo "[step 2] mini_energy (grompp)" | tee -a ./result/b_model.log >&2
gmx mdrun \
    -s ./build/mini_re.tpr \
    -deffnm ./build/mini_re \
    -ntmpi 1 \
    -ntomp "$NPOS" \
    -v >> ./result/b_model.log  2>&1 \
    || { echo -e "${ERROR} gmx mdrun failed (mini_energy)${NC}" | tee -a ./result/b_model.log >&2; exit 1; }

echo -e "${GREEN}Successfully performed a balanced simulation..${NC}" | tee -a ./result/b_model.log >&2
