#!/bin/bash

# 检查是否提供了参数
if [ $# -eq 0 ]; then
    echo "use: $0 [angD|angle]"
    exit 1
fi

# 获取第一个参数作为命令
command="$1"

# 根据输入的命令执行不同操作
case "$command" in
    angD)
        echo "calculate the angle of dipole ang Z axis..."
        sh /data1/huangnan/code/cmd_bash/get_angle.sh
        ;;
    angle)
        echo "calculate angle of 2 vector..."
        sh xxx
        ;;
    *)
        echo "not found command: $command"
        echo "angD|angle can be used"
        exit 1
        ;;
esac