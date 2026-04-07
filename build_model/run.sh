#!/bin/bash

# 当前脚本所在目录，后面调用同目录脚本时要用到
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# 断点续跑标记文件：每完成一步就往这里写一个标记
step_cpt="${bash_dir}/.step.cpt"
# 日志文件
result_log="./result/b_model.log"

# 追加一步完成标记，供下次恢复运行时判断
append_step_marker() {
    local marker="$1"
    echo "$marker" >> "$step_cpt"
}

# 检测应该从哪一步继续：
# 1 表示从头开始
# 2 表示第一步已完成
# 3 表示前两步已完成
# 4 表示全部步骤都已完成
detect_resume_step() {
    local step=1
    if [[ -f "$step_cpt" ]]; then
        if grep -q "^STEP_GET_FRAMES_DONE$" "$step_cpt"; then
            step=4
        elif grep -q "^STEP_EQUILIBRATION_DONE$" "$step_cpt"; then
            step=3
        elif grep -q "^STEP_BUILD_MODEL_DONE$" "$step_cpt"; then
            step=2
        fi
    fi
    echo "$step"
}

# 根据断点文件判断当前是全新运行还是断点恢复
resume_step=$(detect_resume_step)
if [[ "$resume_step" -eq 1 ]]; then
    # 从头开始时，清理旧的中间目录和结果目录
    rm -rf bulk
    rm -rf build
    rm -rf nvt20
    rm -rf initial
    rm -rf result
    mkdir -p result
    # 清空断点文件
    : > "$step_cpt"
else
    # 断点恢复时保留已有中间文件，并提示从哪一步继续
    mkdir -p result
    echo -e "${YELLOW}[INFO]Checkpoint detected in ${step_cpt}, resume from step ${resume_step}.${NC}" | tee -a "$result_log" >&2
fi

# 检查 GROMACS 环境配置文件是否存在，存在才继续执行
if [[ ! -f "$GMXRC" ]]; then
    echo -e "${RED}[ERROR]GMXRC file not found: $GMXRC${NC}" | tee -a "$result_log" >&2
    exit 1
fi
source "$GMXRC"

# 如果断点记录显示全部完成，就直接退出
if [[ "$resume_step" -ge 4 ]]; then
    echo -e "${GREEN}[OK]All build_model steps were already completed according to ${step_cpt}.${NC}" | tee -a "$result_log" >&2
    exit 0
fi

# 如果没有手动提供 set_density，并且当前是从第 1 步开始，
# 就自动调用脚本计算 bulk density
if [ -z "${set_density:-}" ] && [[ "$resume_step" -le 1 ]]; then
    output=$(sh "${bash_dir}/get_bulk_density.sh")

    # 要求脚本执行成功，并且输出必须包含 OUTPUT: 这一行
    if [ $? -ne 0 ] || ! echo "$output" | grep -q "^OUTPUT:"; then
        echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    else
        echo -e "${OK}success get bulk density ${output}" | tee -a "$result_log" >&2
        # 从输出中提取密度值
        density=$(echo "$output" | grep "^OUTPUT:" | cut -d' ' -f2)
        # 校验提取出来的是不是合法数字
        if ! [[ "$density" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
            echo -e "${ERROR}Error: Invalid density value returned: $density${NC}" | tee -a "$result_log" >&2
            echo -e "${ERROR} get_bulk_density.sh failed${NC}" | tee -a "$result_log" >&2
            exit 1
        else
            echo -e "${OK}success WRITE bulk density ${density}" | tee -a "$result_log" >&2
        fi
    fi
    # 统一格式化为 10 位小数，并导出给后续脚本使用
    density=$(printf "%.10f" "$density")
    export set_density="${density}"
elif [ -z "${set_density:-}" ] && [[ "$resume_step" -gt 1 ]]; then
    # 断点恢复时如果 set_density 为空，则默认继续沿用已有中间结果
    echo -e "${YELLOW}[WARNING]set_density is empty, but workflow is resuming from step ${resume_step}. Continue with existing intermediate files.${NC}" | tee -a "$result_log" >&2
fi

# 第一步：构建模型
if [[ "$resume_step" -le 1 ]]; then
    if sh "${bash_dir}/build_model.sh"; then
        append_step_marker "STEP_BUILD_MODEL_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  I. success! build_model.sh completed" >&2
        echo "  The bulk density has been equilibrated." >&2
        echo "##############################################" >&2
        echo -e "${NC}"
        # 删除这一步产生的临时 pdb 文件
        rm -f step*.pdb
    else
        echo -e "${ERROR} build_model.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step I (build_model.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

# 第二步：平衡化
if [[ "$resume_step" -le 2 ]]; then
    if sh "${bash_dir}/equilibration.sh"; then
        append_step_marker "STEP_EQUILIBRATION_DONE"
        echo -e "${GREEN}"
        echo "##############################################" >&2
        echo "  II. success! equilibration.sh completed" >&2
        echo "##############################################" >&2
        echo -e "${NC}"
    else
        echo -e "${ERROR} equilibration.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step II (equilibration.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi

# 第三步：提取轨迹帧并整理最终输出
if [[ "$resume_step" -le 3 ]]; then
    if sh "${bash_dir}/get_frames.sh"; then
        append_step_marker "STEP_GET_FRAMES_DONE"
        # 复制最终拓扑文件到结果文件名
        cp "${TOP}.top" result.top
        # 如果存在备份拓扑文件，则恢复原始文件名
        if [[ -f baktop1.top ]]; then
            mv baktop1.top "${TOP}.top"
        fi
        echo ""
        echo -e "\033[1;32mALL STEPS COMPLETED SUCCESSFULLY\033[0m" >&2
        echo ""
    else
        echo -e "${ERROR} get_frames.sh failed${NC}" | tee -a "$result_log" >&2
        exit 1
    fi
else
    echo -e "${YELLOW}[INFO]Skip Step III (get_frames.sh), already marked done.${NC}" | tee -a "$result_log" >&2
fi
