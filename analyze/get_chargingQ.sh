#!/bin/bash

# 函数：打印帮助信息
print_help() {
    cat << 'EOF'
Usage: ./run.sh [begin_case] [end_case] [save_str] [sel_ener]

Arguments:
  begin_case     起始 case 编号 (default: 1)
  end_case       结束 case 编号 (default: 1)
  save_str       保存的 .mat 文件名 (default: 自动从路径生成)
  sel_ener       能量选择模式：
                   0 -> "Total-Energy" 和 "Conserved-En" (NVT)
                   1 -> 能量项编号列表 "1 2 3 ... 13" (组分解构)

Auto-generated save_str:
  如果未提供 save_str, 将从当前路径提取关键信息生成。
  示例：
    PC/.../2V/.../0ps/  ->  PC_scan2V0ps.mat
    NC/.../5V/.../100ps/ -> NC_scan5V100ps.mat

Examples:
  ./run.sh 1 10 "slow.mat" 0
  ./run.sh                        # 使用所有默认值
  ./run.sh --help                # 显示此帮助信息

Note:
  MATLAB script './matlab/getQ.m' will be executed.
EOF
}

# 检查是否请求帮助
for arg in "$@"; do
    if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then
        print_help
        exit 0
    fi
done

# === mrun 函数：运行 MATLAB 脚本 ===
mrun() {
    fnm=${1%.*}
    /home/apps/matlab2022/bin/matlab -nodesktop -nosplash -r $fnm
}

# === 自动生成 save_str 的逻辑 ===
generate_save_str() {
    local pwd_path
    pwd_path=$(pwd)

    # 分割路径为数组
    IFS='/' read -ra parts <<< "$pwd_path"

    local pc_part=""
    local v_part=""
    local ps_part=""

    for part in "${parts[@]}"; do
        [[ -z "$part" ]] && continue

        # 提取材料类型：PC、NC、Graphene、CNT（只取第一个）
        if [[ -z "$pc_part" && ( "$part" == "PC" || "$part" == "NC" || "$part" == "Graphene" || "$part" == "CNT" ) ]]; then
            pc_part="$part"
        fi

        # 提取电压：2V, 5V
        if [[ -z "$v_part" && "$part" =~ ^[0-9]+V$ ]]; then
            v_part="$part"
        fi

        # 提取时间：0ps, 100ps
        if [[ -z "$ps_part" && "$part" =~ ^[0-9]+ps$ ]]; then
            ps_part="$part"
        fi
    done

    # 构造文件名
    if [[ -n "$pc_part" && -n "$v_part" && -n "$ps_part" ]]; then
        echo "${pc_part}_scan${v_part}${ps_part}.mat"
    elif [[ -n "$pc_part" ]]; then
        echo "${pc_part}_scan_auto.mat"
    else
        echo "scan_auto.mat"
    fi
}

# === 参数解析 ===
begin_case=${1:-1}
end_case=${2:-1}
save_str_input=${3:-}           # 用户输入的 save_str，为空表示没给
sel_ener=${4:-0}                # 默认为 0

# === 导出基本变量 ===
export begin_case
export end_case

# === 生成 save_str：用户没给就自动生成 ===
if [[ -z "$save_str_input" ]]; then
    save_str=$(generate_save_str)
    echo "自动检测路径并生成 save_str: $save_str"
else
    save_str="$save_str_input"
fi
export save_str

# === 设置 selected_string ===
if [[ "$sel_ener" -eq 1 ]]; then
    selected_string="1 2 3 4 5 6 7 8 9 10 11 12 13"
elif [[ "$sel_ener" -eq 0 ]]; then
    selected_string=$'Total-Energy\nConserved-En'  # 使用 $'' 才能解析 \n
else
    echo "警告: sel_ener 应为 0 或 1, 当前值为 $sel_ener, 默认使用 Total-Energy 和 Conserved-En" >&2
    selected_string=$'Total-Energy\nConserved-En'
fi
export selected_string
export save_str
export selected_string

echo "参数配置:"
echo "  begin_case = $begin_case"
echo "  end_case   = $end_case"
echo "  save_str   = $save_str"
echo "  sel_ener   = $selected_string"
echo "running..."

cp /data1/huangnan/code/matlab/getQ.m .
mrun getQ.m