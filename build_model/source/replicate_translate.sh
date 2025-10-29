#!/bin/bash

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if gmx editconf -f ./build/pre_eq.gro -o ./build/pre_eq2.gro -translate 0 0 ${zbox} >> ./result/b_model.log 2>&1;then
    echo -e "${OK}Successfully translated structure along z-axis by ${zbox} nm." | tee -a ./result/b_model.log >&2
else
    echo -e "${ERROR}Translation failed.${NC}" | tee -a ./result/b_model.log >&2
    exit 1
fi

cpp_src=${bash_dir}/merge_gro.cpp
execu_bin=${bash_dir}/../../bin/replicate_translate
bash ${bash_dir}/build_cpp.sh ${cpp_src} ${execu_bin}

${execu_bin} ./build/pre_eq.gro ./build/pre_eq2.gro

double_z_box=$(echo "$zbox * 2" | bc -l)
formatted_line=$(printf "%8.4f%8.4f%8.4f" $xbox $ybox $double_z_box)
sed -i "\$ s/^.*$/      $formatted_line/" ./build/pre_eq_merge.gro
echo -e "${OK}${GREEN}merge two gro successfully!${NC}" | tee -a ./result/b_model.log >&2

cpp_src=${bash_dir}/double_slit_top.cpp
execu_bin=${bash_dir}/../../bin/double_slit_top
bash ${bash_dir}/build_cpp.sh ${cpp_src} ${execu_bin}
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
cp "${MDP}".mdp ./build/grompp_rebalance.mdp
cp "${mini_MDP}".mdp ./build/mini.mdp
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_reB}/" ./build/grompp_rebalance.mdp
sed -i -e 's/^freezegrps.*$/freezegrps = CL  CR  GRA/' -e 's/^freezedim.*$/freezedim   = Y Y Y Y Y Y Y Y Y/' ./build/grompp_rebalance.mdp
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


echo "[step 3] genenrate .tpr  (grompp)" | tee -a ./result/b_model.log >&2
gmx grompp \
    -f ./build/grompp_rebalance.mdp \
    -c ./build/mini_re.gro \
    -o ./build/pre_eq.tpr \
    -maxwarn "${maxWarn}" \
    -p "${TOP}.top" \
    -n ./build/index.ndx  >> ./result/b_model.log  2>&1 \
    || { echo -e "${ERROR} gmx grompp failed${NC}" | tee -a ./result/b_model.log >&2; exit 1; }

echo "[step 4]  mdrun" | tee -a ./result/b_model.log >&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\tUsing GPU acceleration"  | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./build/pre_eq.tpr
        -deffnm ./build/pre_eq
        -ntmpi 1
        -ntomp "$NPOS"
        -pme gpu
        -pmefft gpu
        -nb gpu
        -tunepme no
        -v
    )
else
    echo -e "\tUsing CPU computation"  | tee -a ./result/b_model.log >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./build/pre_eq.tpr
        -deffnm ./build/pre_eq
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi

echo -e "\t${mdrun_cmd[*]}" | tee -a ./result/b_model.log >&2
"${mdrun_cmd[@]}" >> ./result/b_model.log  2>&1
echo -e "${GREEN}Successfully performed a balanced simulation..${NC}" | tee -a ./result/b_model.log >&2
