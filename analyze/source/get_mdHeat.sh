#!/bin/bash
set -euo pipefail

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

g++ -o ./get_mdHeat ${bash_dir}/source/get_mdHeat.cpp -O3 

./get_mdHeat

echo -e "${GREEN}MdHeat data calculation completed!${NC}"