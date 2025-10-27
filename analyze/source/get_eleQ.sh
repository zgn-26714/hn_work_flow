#!/bin/bash
set -euo pipefail

workdir="$(pwd)"
bash_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

mkdir -p deal_data/eleQ/
if g++ -o get_eleQ ${bash_dir}/module_C++/get_eleQ.cpp -O3 ; then #可能需要增加额外的静态链接库
    ./get_eleQ
    echo -e "${GREEN}EleQ data calculation completed!${NC}"
else
    echo -e "${ERROR}complie get_eleQ.cpp failed!${NC}"
    exit 1
fi