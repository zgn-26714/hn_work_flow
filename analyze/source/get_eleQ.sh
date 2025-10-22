#!/bin/bash
set -euo pipefail

workdir="$(pwd)"
bash_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

mkdir -p deal_data
g++ -o get_eleQ ${bash_dir}/module_C++/get_eleQ.cpp -O3
./get_eleQ