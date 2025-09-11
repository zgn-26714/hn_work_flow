#!/bin/bash

mrun() {
    fnm=${1%.*}
    /home/apps/matlab2022/bin/matlab -nodesktop -nosplash -r $fnm
}

# === 显示帮助信息 ===
show_help() {
    cat << 'EOF'
Usage:
    ./script.sh [num] [begin_n] [isD] [save_str] [isGeo] [out_str] [begin_t] [end_t] [onfly_in]

Description:
    Automates the execution of the onflyPost_dens3D_dynamic tool and invokes MATLAB for post-processing.
    All arguments are optional and default values will be used if omitted.

Arguments:
    num             : Number of files (default: 1)
    begin_n         : Starting folder index (default: 1)
    isD             : Dynamic mode flag (0=static, 1=dynamic) (default: 1)
    Environment Variables Passed to MATLAB:
        save_str    : Output .mat file name
        isGeo       : Whether to use geometric center-based number density (0 or 1)
        out_str     : Output file name for onfly (if specified)
    begin_t         : Starting timestep (default: 0)
    end_t           : Ending timestep (default: 0, means unlimited)
    onfly_in        : Name of the *.onfly input file (default: densMNC3D)

Options:
    -h, --help      : Show this help message and exit
    --onflyH        : Print help information from the onflyPost_dens3D_dynamic executable

Examples:
    ./script.sh                             # Use all default values
    ./script.sh 10                          # Set num=10
    ./script.sh 10 1 1 "" 1 output.dat 0 100
                                            # Run with all parameters (save_str empty => auto-generated)
    ./script.sh --help                      # Display help

Notes:
    - If save_str is empty or not provided, it will be auto-generated based on the current path (e.g., PC5nm_scan2V100ps.mat).
    - If out_str is empty, onflyPost_dens3D_dynamic will use its default output filename.
    - Temporary files will be automatically cleaned up after MATLAB script finishes.
EOF
}

# === 自动生成 save_str 的逻辑 ===
generate_save_str() {
    local pwd_path
    pwd_path=$(pwd)

    local pc_part=""
    local v_part=""
    local ps_part=""
    local nm_part=""

    IFS='/' read -ra parts <<< "$pwd_path"

    for part in "${parts[@]}"; do
        [[ -z "$part" ]] && continue

        # 材料类型
        if [[ -z "$pc_part" && ( "$part" == "PC" || "$part" == "ACN" || "$part" == "SOLwork" || "$part" == "DME" ) ]]; then
            pc_part="$part"
        fi

        # 尺寸：5nm, 10nm
        if [[ -z "$nm_part" && "$part" =~ ^[0-9]+nm$ ]]; then
            nm_part="$part"
        fi

        # 电压：2V, 5V
        if [[ -z "$v_part" && "$part" =~ ^[0-9]+V$ ]]; then
            v_part="$part"
        fi

        # 时间：100ps, 200ps
        if [[ -z "$ps_part" && "$part" =~ ^[0-9]+ps$ ]]; then
            ps_part="$part"
        fi
    done

    # 构造文件名
    if [[ -n "$pc_part" && -n "$v_part" && -n "$ps_part" && -n "$nm_part" ]]; then
        echo "${pc_part}${nm_part}_scan${v_part}${ps_part}.mat"
    elif [[ -n "$pc_part" && -n "$v_part" && -n "$ps_part" ]]; then
        echo "${pc_part}_scan${v_part}${ps_part}.mat"
    elif [[ -n "$pc_part" ]]; then
        echo "${pc_part}_scan.mat"
    else
        echo "scan.mat"
    fi
}

# === 参数解析 ===
# 检查帮助选项
for arg in "$@"; do
    case "$arg" in
        -h|--help)
            show_help
            exit 0
            ;;
        --onflyH)
            onflyPost_dens3D_dynamic -h
            exit 0
            ;;
    esac
done

# 设置默认值
num=${1:-1}
begin_n=${2:-1}
isD=${3:-1}
save_str_input=${4:-}     # 用户指定的 save_str（MATLAB 保存文件）
isGeo=${5:-0}             # 几何参数，传给 MATLAB
out_str=${6:-}            # onfly 输出文件名
begin_t=${7:-0}
end_t=${8:-0}
onfly_in=${9:-"densMNC3D"}

# === 自动生成 save_str（如果未提供）===
if [[ -z "$save_str_input" ]]; then
    save_str=$(generate_save_str)
    echo "=> 自动生成 save_str: $save_str"
else
    save_str="$save_str_input"
fi

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

# === 导出环境变量供 MATLAB 使用 ===
export save_str
export isGeo

# === 运行 onflyPost_dens3D_dynamic ===
echo "=> onflyPost_dens3D_dynamic..."
if [[ -z "$out_str" ]]; then
    echo "commmand: onflyPost_dens3D_dynamic "../case$begin_n/$onfly_in.onfly" -n "$num" -d "$isD" -b "$begin_t" -e "$end_t" "
    export out_str="onfly3D.dat"
    onflyPost_dens3D_dynamic "./case$begin_n/$onfly_in.onfly" -n "$num" -d "$isD" -b "$begin_t" -e "$end_t"
else
    echo "commmand: onflyPost_dens3D_dynamic "../case$begin_n/$onfly_in.onfly" -n "$num" -d "$isD" -o "$out_str" -b "$begin_t" -e "$end_t" "
    export out_str
    onflyPost_dens3D_dynamic "./case$begin_n/$onfly_in.onfly" -n "$num" -d "$isD" -o "$out_str" -b "$begin_t" -e "$end_t"
fi

if [ $? -ne 0 ]; then
    echo "=> 错误: onflyPost_dens3D_dynamic 执行失败！"
    exit 1
fi

echo "running remake_onfly.m..."
# === 复制 MATLAB 脚本 ===
cp "/data1/huangnan/code/matlab/remake_onfly.m" ./
if [ $? -ne 0 ]; then
    echo "=> 错误: 无法复制 remake_onfly.m 文件！"
    exit 1
fi

# === 运行 MATLAB 脚本 ===
mrun "remake_onfly.m"

# === 清理临时文件 ===
rm -f "remake_onfly.m"
echo "=> 临时文件 remake_onfly.m 已删除。"

echo "=================================================="
echo "           所有任务执行完毕！"
echo "=================================================="