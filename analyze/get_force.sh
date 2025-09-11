#!/bin/bash

# 设置循环参数
begin=$1
last=$2

echo "开始执行test_force循环，从case${begin}到case${last}"
echo "=========================================="

# 循环执行命令
for ((i=$begin; i<=$last; i++)); do
    echo "正在处理 case${i}..."
    
    # 检查必要文件是否存在
    case_dir="./case${i}"
    xtc_file="${case_dir}/slownvt.xtc"
    tpr_file="${case_dir}/slownvt.tpr"
    ndx_file="${case_dir}/index.ndx"
    
    if [[ ! -d "$case_dir" ]]; then
        echo "错误: 目录 $case_dir 不存在，跳过..."
        continue
    fi
    
    if [[ ! -f "$xtc_file" ]]; then
        echo "警告: 文件 $xtc_file 不存在"
    fi
    
    if [[ ! -f "$tpr_file" ]]; then
        echo "警告: 文件 $tpr_file 不存在"
    fi
    
    if [[ ! -f "$ndx_file" ]]; then
        echo "警告: 文件 $ndx_file 不存在"
    fi
    
    # 执行test_force命令
    echo "执行命令: test_force -f $xtc_file -s $tpr_file -n $ndx_file -n2ave 100 -o ./xvgdat/cal_bond_force/"$i"force21_29.xvg -oset 0 -zmin 21 -zmax 29"
    
    test_force -f "$xtc_file" -s "$tpr_file" -n "$ndx_file" -n2ave 100 -o ./xvgdat/cal_bond_force/"$i"force21_29.xvg -oset 0 -zmin 21 -zmax 29
    
    # 检查命令执行状态
    if [[ $? -eq 0 ]]; then
        echo "✓ case${i} 执行成功"
    else
        echo "✗ case${i} 执行失败 (退出码: $?)"
    fi
    
    echo "----------------------------------------"
done

echo "所有任务完成！"