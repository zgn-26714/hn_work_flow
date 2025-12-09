#!/bin/bash
set -euo pipefail

outfile="${packmol}.inp"   # 输出文件

read -r -a MOL_arr <<< "$MOL"
read -r -a MOLnum_arr <<< "$MOL_num"

last_line=$(tail -n 1 ./model/single_electrode.gro)

read -r BOX_X BOX_Y BOX_Z <<< "$last_line"



cat >> "$outfile" << EOF
structure single_electrode.pdb
number 1
fixed  0. 0. 0. 0. 0. 0.
end structure
EOF

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
