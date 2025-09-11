#!/bin/bash
echo "">&2
echo "**********************************************************************">&2
echo "     3/3   get frames">&2
echo "⚙️  Parameters: nsteps = ${nsteps_ini}">&2
echo "**********************************************************************">&2
if [[ -f "$GMXRC" ]]; then
    source "$GMXRC"
else
    echo "【ERROR】GMXRC not found: $GMXRC" >&2
    exit 1
fi

echo "【step 1】 clear file"
rm -f *.tpr
rm -f *#
rm -rf initial
mkdir -p initial
echo "【SUCCESS】">&2

sed -i "s/^[[:space:]]*nsteps[[:space:]]*=.*/nsteps = ${nsteps_ini}/" "$MDP".mdp
old_nstxtcout=$(grep -E '^[[:space:]]*nstxtcout[[:space:]]*=' "$MDP".mdp | sed -E 's/.*=[[:space:]]*([0-9]+).*/\1/')
sed -i -E \
    -e "s/^([[:space:]]*nstxtcout[[:space:]]*=[[:space:]]*)[0-9]+/\10/" \
    -e "s/^([[:space:]]*nstxout[[:space:]]*=[[:space:]]*)[0-9]+/\1$old_nstxtcout/" \
    "$MDP".mdp

echo "【step 2】 genenrate .tpr (grompp)">&2
gmx grompp \
    -f "$MDP".mdp \
    -c ./nvt20/nvt20.gro \
    -o ./initial/initial.tpr \
    -maxwarn "${maxWarn}" \
    -p "${TOP}.top" \
    -n "${INDEX}.ndx" > test.log 2>&1 \
    || { echo "【ERROR】gmx grompp failed" >&2; exit 1; }
echo "【SUCCESS】">&2

#  mdrun 
echo "【step 3】  mdrun">&2
if [[ "$GPU" -eq 1 ]]; then
    echo -e "\t【INFO】Using GPU acceleration" >&2
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
    echo -e "\t【INFO】Using CPU computation" >&2
    mdrun_cmd=(
        gmx mdrun
        -s ./initial/initial.tpr
        -deffnm ./initial/initial
        -ntmpi 1
        -ntomp "$NPOS"
        -v
    )
fi
echo -e "\t【CMD】${mdrun_cmd[*]}">&2
"${mdrun_cmd[@]}" > test.log 2>&1 
echo "【SUCCESS】">&2

echo "【step 4】 get frames" >&2

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
get_frame_bin="${bash_dir}/source/get_frame"

# 检查可执行文件是否存在
if [[ ! -x "${get_frame_bin}" ]]; then
    echo "【INFO】'get_frame' not found or not executable. Attempting to compile..." >&2
    src_file="${bash_dir}/source/get_frame.cpp"  # 假设源文件是 get_frame.cpp
    if [[ ! -f "${src_file}" ]]; then
        echo "【ERROR】Source file '${src_file}' not found. Cannot compile." >&2
        exit 1
    fi

    # 尝试编译
    g++ -O3 -std=c++17 -o "${get_frame_bin}" "${src_file}"
    if [[ $? -ne 0 ]]; then
        echo "【ERROR】Compilation of 'get_frame' failed." >&2
        exit 1
    fi

    # 确保编译后文件可执行
    chmod +x "${get_frame_bin}"
    echo "【INFO】'get_frame' compiled successfully." >&2
fi

# 执行程序
"${get_frame_bin}" -i ./initial/initial -t "${Temperature}" -n "${num_intialframes}" -o "${Temperature}k"
if [[ $? -ne 0 ]]; then
    echo "【ERROR】Execution of 'get_frame' failed." >&2
    exit 1
fi

echo "【SUCCESS】" >&2