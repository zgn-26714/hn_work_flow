#!/bin/bash
rm ./result/deal_data.log

set -euo pipefail
# 检查是否提供了参数
if [ $# -eq 0 ]; then
    echo -e "${ERROR} Missing parameter for analyze command.${NC}" | tee -a ./result/deal_data.log >&2
    echo "use: $0 [JE|...]"  | tee -a ./result/deal_data.log >&2
    exit 1
fi

mrun(){
    fnm=${1%.*}
    /home/apps/matlab2022/bin/matlab -nodesktop -nosplash -r $fnm
}

export -f mrun


command="$1"
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# 根据输入的命令执行不同操作
case "$command" in
    JE)
        export onfly_dat=$2
        echo "Calculate and store the data required for the JE model..."  | tee -a ./result/deal_data.log >&2
        bash ${bash_dir}/source/JE_model.sh  | tee -a ./result/deal_data.log >&2
        ;;

    *)
        echo -e "${ERROR}not found command: $command" | tee -a ./result/deal_data.log >&2
        exit 1
        ;;
esac