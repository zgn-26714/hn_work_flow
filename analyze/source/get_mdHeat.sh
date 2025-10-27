#!/bin/bash
set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

mkdir -p ./deal_data/mdHeat/
#g++ -o ./get_mdHeat ${bash_dir}/module_C++/get_mdHeat.cpp -O3 ${bash_dir}/module_C++/include/libmatio.a -I${bash_dir}/module_C++/include/ -lm -lz -lpthread
g++ -o ./get_mdHeat ${bash_dir}/module_C++/get_mdHeat.cpp -O3 \
    ${bash_dir}/module_C++/include/libmatio.a \
    /data1/huangnan/github/matlab/hd5/hdf5-1.14.3/build/hn_install/hdf5/lib/libhdf5.a \
    /data1/huangnan/github/matlab/hd5/hdf5-1.14.3/build/hn_install/hdf5/lib/libhdf5_hl.a \
    /data1/huangnan/github/matlab/hd5/hdf5-1.14.3/build/hn_install/hdf5/lib/libhdf5_tools.a \
    -I${bash_dir}/module_C++/include/ \
    -I/data1/huangnan/github/matlab/hd5/hdf5-1.14.3/build/hn_install/hdf5/include \
    -lz -lm -lpthread -ldl
./get_mdHeat

echo -e "${GREEN}MdHeat data calculation completed!${NC}"