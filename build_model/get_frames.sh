#!/bin/bash

{
echo -e "${BLUE}"
echo "**********************************************************************">&2
echo "     3/3   get frames">&2
echo "⚙️  Parameters: nsteps = ${nsteps_ini}">&2
echo "**********************************************************************">&2
echo -e "${NC}"
} | tee -a ./result/b_model.log 
if [[ -f "$GMXRC" ]]; then
    source "$GMXRC"
else
    echo -e "${ERROR}GMXRC not found: $GMXRC${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi

echo "[step 1] clear file" | tee -a ./result/b_model.log
rm -f *#

mkdir -p initial
echo -e "${OK}" | tee -a ./result/b_model.log >&2

cp "${MDP}.mdp" ./initial/grompp.mdp
sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_ini}/" ./initial/grompp.mdp

sed -i -E \
    -e "s/^([[:space:]]*nstxtcout[[:space:]]*=[[:space:]]*)[0-9]+/\10/" \
    -e "s/^([[:space:]]*nstxout[[:space:]]*=[[:space:]]*)[0-9]+/\1${XOUT_FRAMES}/" \
    ./initial/grompp.mdp

echo "[step 2] genenrate .tpr (grompp)" | tee -a ./result/b_model.log >&2
echo q | gmx make_ndx -f ./nvt20/nvt20.gro -o ./initial/index.ndx >> ./result/b_model.log  2>&1

gmx grompp \
    -f ./initial/grompp.mdp \
    -c ./nvt20/nvt20.gro \
    -o ./initial/initial.tpr \
    -maxwarn "${maxWarn}" \
    -p "${TOP}.top" \
    -n ./initial/index.ndx >> ./result/b_model.log  2>&1 \
    || { echo -e "${ERROR}gmx grompp failed${NC}" | tee -a ./result/b_model.log >&2; exit 1; }
echo -e "${OK}" | tee -a ./result/b_model.log >&2

#  mdrun 
echo "[step 3]  mdrun" | tee -a ./result/b_model.log >&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\tUsing GPU acceleration" | tee -a ./result/b_model.log  >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./initial/initial.tpr
        -deffnm ./initial/initial
        -ntmpi 1
        -ntomp "$NPOS"
        -pme gpu
        -pmefft gpu
        -nb gpu
        -tunepme no
        -v
    )
else
    echo -e "\tUsing CPU computation" | tee -a ./result/b_model.log  >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./initial/initial.tpr
        -deffnm ./initial/initial
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi
echo -e "\t${mdrun_cmd[*]}" | tee -a ./result/b_model.log >&2
"${mdrun_cmd[@]}" >> ./result/b_model.log  2>&1 
echo -e "${OK}" | tee -a ./result/b_model.log >&2

echo "[step 4] get frames"  | tee -a ./result/b_model.log >&2

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
get_frame_bin=./get_frame

# 检查可执行文件是否存在
if [[ ! -x "${get_frame_bin}" ]]; then
    echo "'get_frame' not found or not executable. Attempting to compile..."  | tee -a ./result/b_model.log >&2
    src_file="${bash_dir}/source/get_frame.cpp"  # 源文件是 get_frame.cpp
    if [[ ! -f "${src_file}" ]]; then
        echo -e "${ERROR}Source file '${src_file}' not found. Cannot compile.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    # 尝试编译
    g++ -O3 -std=c++17 -o "${get_frame_bin}" "${src_file}"
    if [[ $? -ne 0 ]]; then
        echo -e "${ERROR}Compilation of 'get_frame' failed.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    # 确保编译后文件可执行
    chmod +x "${get_frame_bin}"
    echo "'get_frame' compiled successfully." | tee -a ./result/b_model.log  >&2
fi

# 执行程序
echo "${num_intialframes}" >&2
Temperature=$(grep "ref_t" ./initial/grompp.mdp | awk -F'=' '{print $2}' | awk '{print $1}')
"${get_frame_bin}" -i ./initial/initial -t "${Temperature}" -n "${num_intialframes}" -o "${Temperature}k"
if [[ $? -ne 0 ]]; then
    echo -e "${ERROR}Execution of 'get_frame' failed.${NC}" | tee -a ./result/b_model.log  >&2
    exit 1
fi

echo -e "${GREEN}Successfully generated the required initial structure.${NC}" | tee -a ./result/b_model.log >&2