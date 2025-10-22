#!/bin/bash
set -euo pipefail

echo "running remake_onfly.m..."

mrun(){
    fnm=${1%.*}
    /home/apps/matlab2022/bin/matlab -nodesktop -nosplash -r $fnm
}



bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cp "${bash_dir}/module_C++/remake_onfly_den.m" ./
if [ $? -ne 0 ]; then
    echo "=> 错误: 无法复制 remake_onfly.m 文件！"
    exit 1
fi


mrun "remake_onfly_den.m"
rm -f "remake_onfly_den.m"

echo -e "${GREEN}Successfully completed remake_onfly.m!${NC}"
