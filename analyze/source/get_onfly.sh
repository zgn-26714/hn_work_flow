#!/bin/bash
set -euo pipefail

# === 输出参数摘要 ===
echo "=================================================="
echo "       onflyPost_dens3D_dynamic 脚本启动"
echo "--------------------------------------------------"
printf "%-20s : %s\n" "num (文件数量)" "$num"
printf "%-20s : %s\n" "begin_n (起始文件夹索引)" "$begin_n"
printf "%-20s : %s\n" "isD (动态模式)" "$isD"
printf "%-20s : %s\n" "save_str (MATLAB输出)" "$save_str"
printf "%-20s : %s\n" "isGeo (geo_numD)" "$isGeo"
printf "%-20s : %s\n" "out_str (onfly输出)" "${out_str:-default}"
printf "%-20s : %s\n" "begin_t (起始时间步)" "$begin_t"
printf "%-20s : %s\n" "end_t (结束时间步)" "$end_t"
echo "=================================================="


# === 运行 onflyPost_dens3D_dynamic ===
echo "=> onflyPost_dens3D_dynamic..."
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

g++ -o onflyPost_dens3D_dynamic ${bash_dir}/module_C++/onflyPost_densANDBALJ3D_dynamic.cpp -O3

analyze_num=((${analyze_end_case}-${analyze_begin_case}+1))
./onflyPost_dens3D_dynamic "./case${analyze_begin_case}/${analysis_ONFLY_in}.onfly" \
-n "${analyze_num}" -d "${analysis_ONFLY_isD}" -o ./deal_data/onfly/onfly"${analyze_begin_case}-${analyze_end_case}".dat \
-b "${analysis_begin_t}" -e "${analysis_end_t}"


if [ $? -ne 0 ]; then
    echo "=> 错误: onflyPost_dens3D_dynamic 执行失败！"
    exit 1
else
    echo -e "${GREEN}Successfully completed onflyPost_dens3D_dynamic!${NC}"
fi
