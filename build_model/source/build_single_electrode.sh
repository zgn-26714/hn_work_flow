#!/bin/bash
set -euo pipefail


mkdir -p ./model
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cpp_src=${bash_dir}/generate_slit.cpp
execu_bin=${bash_dir}/../../bin/generate_slit

if [[ ! -x "${execu_bin}" ]]; then
    echo "Executable file not found or not executable. Attempting to compile..."  | tee -a ./result/b_model.log >&2
    if [[ ! -f "${cpp_src}" ]]; then
        echo -e "${ERROR}Source file '${cpp_src}' not found. Cannot compile.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    g++ "$cpp_src" -o "$execu_bin" -O3 -std=c++17
    if [[ $? -ne 0 ]]; then
        echo -e "${ERROR}Compilation of Executable file failed.${NC}" | tee -a ./result/b_model.log  >&2
        exit 1
    fi

    chmod +x "${execu_bin}"
    echo "Executable file compiled successfully." | tee -a ./result/b_model.log  >&2
fi


export BASH_DIR=$bash_dir
${execu_bin} ./model/single_electrode.gro

n=$(wc -l < ./model/single_electrode.gro)
echo "n=$n" | tee -a ./result/b_model.log >&2

read min max < <(
awk -v n="$n" 'NR>2 && NR<n {
    if(!i){min=max=$NF;i=1}
    if($NF<min)min=$NF
    if($NF>max)max=$NF
} END{print min, max}' ./model/single_electrode.gro
)

# 记录 min 和 max 到日志
echo "min=$min max=$max" | tee -a ./result/b_model.log >&2

##更改盒子大小和位置
last_line=$(tail -n 1 ./model/single_electrode.gro)
read -r BOX_X BOX_Y BOX_Z <<< "$last_line"
bulk_half=$(echo "scale=3; ${len_bulk} / 2" | bc)
gmx editconf -f ./model/single_electrode.gro -o ./model/single_electrode.gro -translate 0 0 -${min} >> ./result/b_model.log  2>&1
gmx editconf -f ./model/single_electrode.gro -o ./model/single_electrode.gro -translate 0 0 "${bulk_half}" >> ./result/b_model.log  2>&1
min_ele=${bulk_half}
max_ele=$(awk "BEGIN{print $max - $min + $bulk_half}")
now_BOX_Z=$(awk "BEGIN{print $max_ele + $bulk_half}")
formatted_line=$(printf "%8.4f%8.4f%8.4f" $BOX_X $BOX_Y $now_BOX_Z)
sed -i "\$ s/^.*$/      $formatted_line/" ./model/single_electrode.gro
gmx editconf -f ./model/single_electrode.gro -o ./model/single_electrode.pdb >> ./result/b_model.log  2>&1
echo "Box size updated to: $BOX_X x $BOX_Y x $now_BOX_Z nm" | tee -a ./result/b_model.log  >&2
echo -e "${OK} get initial electrode, Z coordinate range: min=${min_ele}, max=${max_ele}${NC}" | tee -a ./result/b_model.log  >&2

export min_ele
export max_ele
#自动生成packmol
if [ "$is_auto_top" -eq 1 ]; then
    if bash ${bash_dir}/auto_packmol.sh >> ./result/b_model.log  2>&1; then
        echo -e "${OK}Successfully generated ${packmol}.inp for single electrode model." | tee -a ./result/b_model.log >&2
    else
        echo -e "${ERROR}Failed to generate ${packmol}.inp for single electrode model.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi
fi


#自动生成初始拓扑
if [ "$is_auto_top" -eq 1 ]; then
    if bash ${bash_dir}/auto_topology.sh >> ./result/b_model.log  2>&1; then
        echo -e "${OK}Successfully generated topology for single electrode model." | tee -a ./result/b_model.log >&2
    else
        echo -e "${ERROR}Failed to generate topology for single electrode model.${NC}" | tee -a ./result/b_model.log >&2
        exit 1
    fi
fi