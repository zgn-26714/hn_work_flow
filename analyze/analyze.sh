#!/bin/bash
set -euo pipefail

# 检查是否提供了参数
if [ $# -eq 0 ]; then
    echo "${ERROR} Missing parameter for analyze command.${NC}" >&2
    echo "use: $0 [angD|angle|...]" >&2
    exit 1
fi

# 获取第一个参数作为命令
command="$1"
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# 根据输入的命令执行不同操作
case "$command" in
    elecharge)
        echo "calculate the angle of dipole ang Z axis..."
        bash ${bash_dir}/source/get_eleQ.sh
        ;;
    onfly)
        echo "calculate onfly data"
        if bash ${bash_dir}/source/get_onfly.sh; then
            bash ${bash_dir}/source/remake_onfly.sh
        else
            echo -e "${ERROR}onfly data calculation failed!${NC}"
            exit 1
        fi
        ;;
    mdheat)
        echo "calculate mdheat data"
        bash ${bash_dir}/source/get_mdHeat.sh
        ;;
    angle)
        echo "calculate angle of 2 vector..."
        bash xxx
        ;;

    angD)
        echo "calculate elecharge data"
        bash xxx
        ;;
    *)
        echo "not found command: $command"
        echo "angD|angle can be used"
        exit 1
        ;;
esac