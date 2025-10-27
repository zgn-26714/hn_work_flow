#!/bin/bash
set -euo pipefail

echo "Starting script execution..."

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
echo "Script directory: ${bash_dir}"

# 复制文件
echo "Copying remake_onfly_den.m file..."
cp "${bash_dir}/source/remake_onfly_den.m" ./

if [ ! -f "remake_onfly_den.m" ]; then
    echo -e "${RED}Error: remake_onfly_den.m file not found after copying!${NC}"
    exit 1
fi

echo "Executing MATLAB script..."
# 运行MATLAB脚本
mrun remake_onfly_den.m
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to execute remake_onfly_den.m!${NC}"
    exit 1
fi

rm -f remake_onfly_den.m
echo "Temporary files removed."

echo -e "${GREEN}Successfully completed remake_onfly_den.m execution!${NC}"
echo "Script finished."