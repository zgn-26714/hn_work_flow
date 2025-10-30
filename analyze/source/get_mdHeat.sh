#!/bin/bash
set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

mkdir -p ./deal_data/mdHeat/
build_dir=${bash_dir}/../../bin/build_cpp.sh
src_cpp=${bash_dir}/module_C++/get_mdHeat.cpp
execu_bin=${bash_dir}/../../bin/get_mdHeat
#g++ -o ./get_mdHeat ${bash_dir}/module_C++/get_mdHeat.cpp -O3 ${bash_dir}/module_C++/include/libmatio.a -I${bash_dir}/module_C++/include/ -lm -lz -lpthread
bash ${build_dir} ${src_cpp} ${execu_bin}

${execu_bin} >> ./result/analyze.log  2>&1

echo -e "${GREEN}MdHeat data calculation completed!${NC}"