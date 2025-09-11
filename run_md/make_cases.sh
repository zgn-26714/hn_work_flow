#!/bin/bash

# 简易计算器函数（用于浮点运算，虽然本脚本未使用浮点）
calc() {
    python -c "print($*)"
}

show_help() {
    cat << 'EOF'
    you need to modify the slowdir before you use it
Usage: ./script.sh START END [NUM] [QUEUE] [MAKECASE] [MAKEJOB] [mode]

Automatically create and submit charging simulation jobs.

Arguments:
  START        Start case number (e.g., 10)
  END          End case number (exclusive, e.g., 20)
  NUM          Number of cases per job [default: 1]
  QUEUE        PBS queue: short / new / fast [default: new]
  MAKECASE     1 to create case dirs, 0 to skip [default: 1]
  MAKEJOB      1 to generate & submit job, 0 to skip [default: 1]
  mode         append to run md with append mode(default: default)

Examples:
  ./script.sh 10 20
  ./script.sh 100 150 3 fast 1 1
  ./script.sh 1 10 1 new 0 1   # only submit job, skip making case

Note:
  - This script uses PBS (qsub) for job submission.
  - Ensure all paths (workdir, framedir, etc.) are correct.
  - On macOS, use gsed instead of sed if needed.
EOF
}

# 参数检查 + help
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

if [[ -z "$1" || -z "$2" ]]; then
    echo "❌ Error: Missing required arguments."
    echo "Use -h or --help for usage."
    exit 1
fi
# ============== 配置参数 ===============
workdir="$(pwd)"
framedir="$workdir/initialframes"
matrixdir="$workdir/matrixdata"
slowdir="/data1/huangnan/scanratepotential"

[[ -d "$framedir" ]] || { echo "❌ Frame dir not found: $framedir"; exit 1; }
[[ -d "$matrixdir" ]] || { echo "❌ Matrix dir not found: $matrixdir"; exit 1; }
[[ -d "$slowdir" ]] || { echo "❌ Slow dir not found: $slowdir"; exit 1; }
[[ -d "$workdir/basicfile" ]] || { echo "❌ Basic file dir missing: $workdir/basicfile"; exit 1; }
[[ -f "$workdir/runbasic.job2" ]] || { echo "❌ Job template missing: $workdir/runbasic.job2"; exit 1; }

V=2                   # 电压值
ic=0                  # 扫描速率（scan rate）
num=${3:-1}                # 每次处理多少个 case
quene=${4:-"new"}           # PBS 队列名
ismakecase=${5:-1}
ismakejob=${6:-1}
mdoe=${7:-"default"}

cas=$1                # 起始 case 编号
endcas=$2             # 结束 case 编号

# =============== 函数定义 ===============
make_case() {
    local ca=$1
    local objdir="$workdir/charging/298k/${V}V/${ic}ps/case$ca"

    mkdir -p "$objdir"
    cp "$framedir/frame$ca.gro" "$objdir/"
    cp "$workdir/basicfile"/* "$objdir/"
    cp "$matrixdir/CPM_ControlFile.dat_0V" "$objdir/CPM_ControlFile.dat"

    # 删除旧文件并创建符号链接
    rm -f "$objdir/allMatrixA.bin" "$objdir/Dphis_control.dat"
    ln -s "$matrixdir/allMatrixA.bin" "$objdir/allMatrixA.bin"
    ln -s "$slowdir/flexible_Dphis_control.dat${V}_${ic}" "$objdir/Dphis_control.dat"

    echo "Created case: $ca at $objdir" >&2

    echo "$objdir"
}

make_job() {
    local start_case=$1
    local objjob="$workdir/jobfile/j${V}V${ic}pscase${start_case}.job"

    mkdir -p "$workdir/jobfile"
    SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
    cp "$SCRIPT_PATH/runbasic.job2" "$objjob"

    # 批量替换 PBS 脚本中的内容
    work_dir=$2

    sed -i " 
        2s/.*/#PBS -q $quene/;
        3s/.*/#PBS -N s_charging${V}V${ic}case${start_case}/;
        /V=/s/.*/V=$V/;
        /j=/s/.*/j=$ic/;
        /k=/s/.*/k=$start_case/;
        /h=/s/.*/h=$num/;
        /workdir=/s/.*/workdir=$work_dir/;
        /mode=/s/.*/mdoe=$mode/
    " "$objjob"

    # 根据队列设置 walltime 和 ppn
    case "$quene" in
        short)  sed -i '6s/.*/#PBS -l walltime=36:00:00/' "$objjob" ;;
        new)    sed -i '5s/.*/#PBS -l nodes=1:ppn=32/' "$objjob" ;;
        fast)   sed -i '5s/.*/#PBS -l nodes=1:ppn=24/' "$objjob" ;;
    esac

    # 提交作业
    qsub "$objjob"
    echo "Submitted job: $objjob"
}

# =============== 主循环 ===============
while (( cas < endcas )); do
    # 创建多个 case
    if (( ismakecase == 1 )); then
        for (( i = 0; i < num; i++ )); do
            objdir=$(make_case $(( cas + i )))
        done
    fi
    rundir="$workdir/charging/298k/${V}V/${ic}ps"
    # 生成并提交作业
    if (( ismakejob == 1 )); then
        make_job "$cas" "$rundir"
    fi

    # 更新 cas
    (( cas += num ))
done