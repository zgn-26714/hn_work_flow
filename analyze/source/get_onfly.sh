#!/bin/bash
set -euo pipefail

mkdir -p ./deal_data/onfly/
analyze_num=$((${analyze_end_case} - ${analyze_begin_case} + 1))
outdir="./deal_data/onfly/onfly"${analyze_begin_case}-${analyze_end_case}".dat"
# === 输出参数摘要 ===
echo "=================================================="
echo "       onflyPost_dens3D_dynamic 脚本启动"
echo "--------------------------------------------------"
printf "%-20s : %s\n" "num (文件数量)" "$analyze_num"
printf "%-20s : %s\n" "begin_n (起始文件夹索引)" ./case"$analyze_begin_case"
printf "%-20s : %s\n" "isD (动态模式)" "$analysis_ONFLY_isD"
printf "%-20s : %s\n" "out_str (onfly输出)" "${outdir}"
printf "%-20s : %s\n" "begin_t (起始时间步)" "$analysis_begin_t"
printf "%-20s : %s\n" "end_t (结束时间步)" "$analysis_end_t"
echo "=================================================="

# === 运行 onflyPost_dens3D_dynamic ===
echo "=> onflyPost_dens3D_dynamic..."
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

onfly_version="onflyPost_denANDBALJ3D_dynamic"

g++ -o ${onfly_version} ${bash_dir}/module_C++/${onfly_version}.cpp -O3

command="./${onfly_version} "./case${analyze_begin_case}/${analysis_ONFLY_in}.onfly" \
-n "${analyze_num}" -d "${analysis_ONFLY_isD}" -o ./deal_data/onfly/onfly"${analyze_begin_case}-${analyze_end_case}".dat \
-b "${analysis_begin_t}" -e "${analysis_end_t}" -om "${MODE_ONFLY}""


if eval $command; then
    echo -e " ${EEROR} onflyPost_dens3D_dynamic fialed!"
    exit 1
else
    echo -e "${GREEN}Successfully completed onflyPost_dens3D_dynamic!${NC}"
fi
