#!/bin/bash
rm ./result/analyze.log

set -euo pipefail
# 检查是否提供了参数
if [ $# -eq 0 ]; then
    echo -e "${ERROR} Missing parameter for analyze command.${NC}" | tee -a ./result/analyze.log >&2
    echo "use: $0 [eleQ|onfly|mdheat|MD_PP|...]"  | tee -a ./result/analyze.log >&2
    exit 1
fi

# 获取第一个参数作为命令
command="$1"
bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# 根据输入的命令执行不同操作
case "$command" in
    eleQ)
        echo "Calculating CPM electrode charge..."  | tee -a ./result/analyze.log >&2
        bash ${bash_dir}/source/get_eleQ.sh  | tee -a ./result/analyze.log >&2
        ;;
    onfly)
        echo "calculate onfly data..."  | tee -a ./result/analyze.log >&2
        bash ${bash_dir}/source/get_onfly.sh  | tee -a ./result/analyze.log >&2
        ;;
    mdheat)
        echo "calculate mdheat data..."  | tee -a ./result/analyze.log >&2
        bash ${bash_dir}/source/get_mdHeat.sh  | tee -a ./result/analyze.log >&2
        ;;
    MD_PP)
        echo "Run the MD post-processing module..." | tee -a ./result/analyze.log >&2
        bash ${bash_dir}/source/MD_PP.sh | tee -a ./result/analyze.log >&2
        ;;
    *)
        echo -e "${ERROR}not found command: $command" | tee -a ./result/analyze.log >&2
        #echo "angD|angle can be used" | tee -a ./result/analyze.log >&2
        exit 1
        ;;
esac