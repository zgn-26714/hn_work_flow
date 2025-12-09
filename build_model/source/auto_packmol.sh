#!/bin/bash
set -euo pipefail

outfile="${packmol}.inp"   # 输出文件

read -r -a MOL_arr <<< "$MOL"
read -r -a MOLnum_arr <<< "$MOL_num"

last_line=$(tail -n 1 ./model/single_electrode.gro)

read -r BOX_X BOX_Y BOX_Z <<< "$last_line"

echo "MOL is ${MOL_arr[*]}" | tee -a ./result/b_model.log >&2 # 输出 MOL_arr 的值
echo "MOLnum is ${MOLnum_arr[*]}" | tee -a ./result/b_model.log >&2 # 输出 MOL_arr 的值
echo "BOX dimensions are: $BOX_X, $BOX_Y, $BOX_Z" | tee -a ./result/b_model.log >&2 # 输出 MOL_arr 的值


max_ele=$(printf "%.3f" "$(bc -l <<< "$max_ele * 10")")
BOX_X=$(printf "%.3f" "$(bc -l <<< "$BOX_X * 10")")
BOX_Y=$(printf "%.3f" "$(bc -l <<< "$BOX_Y * 10")")
BOX_Z=$(printf "%.3f" "$(bc -l <<< "$BOX_Z * 10")")
min_ele=$(printf "%.3f" "$(bc -l <<< "$min_ele * 10")")
bulk_xbox=$(printf "%.3f" "$(bc -l <<< "$bulk_xbox * 10")")
bulk_ybox=$(printf "%.3f" "$(bc -l <<< "$bulk_ybox * 10")")
bulk_zbox=$(printf "%.3f" "$(bc -l <<< "$bulk_zbox * 10")")


cat >> "$outfile" << EOF
structure single_electrode.pdb
number 1
fixed  0. 0. 0. 0. 0. 0.
end structure

EOF

echo -e "${OK}Added single_electrode structure to Packmol input." | tee -a ./result/b_model.log >&2
# ---------------------------
# 上半部分（z = max_ele）
# ---------------------------
for i in "${!MOL_arr[@]}"; do
    cat >> "$outfile" << EOF
structure ${MOL_arr[$i]}.pdb
number ${MOLnum_arr[$i]}
inside box 0. 0. $max_ele  $BOX_X  $BOX_Y  $BOX_Z
end structure

EOF
done

echo -e "${OK} Added upper region structures to Packmol input." | tee -a ./result/b_model.log >&2

# ---------------------------
# 下半部分（z = 0 → min_ele）
# ---------------------------
for i in "${!MOL_arr[@]}"; do
    cat >> "$outfile" << EOF
structure ${MOL_arr[$i]}.pdb
number ${MOLnum_arr[$i]}
inside box 0. 0. 0.  $BOX_X  $BOX_Y  $min_ele
end structure

EOF
done

echo -e "${OK}Added lower region structures to Packmol input." | tee -a ./result/b_model.log >&2

#######生成bulk_packmol.inp
outfile="${bulk_pa}.inp"   # 输出文件
read -r -a MOLnum_arr_bulk <<< "$MOL_num_bulk"


for i in "${!MOL_arr[@]}"; do
    cat >> "$outfile" << EOF
structure ${MOL_arr[$i]}.pdb
number ${MOLnum_arr_bulk[$i]}
inside box 0. 0. 0. $bulk_xbox  $bulk_ybox  $bulk_zbox
end structure

EOF
done

echo -e "${OK}Generated ${bulk_pa}.inp for bulk system." | tee -a ./result/b_model.log >&2