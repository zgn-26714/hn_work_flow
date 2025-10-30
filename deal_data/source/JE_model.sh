#!/bin/bash
set -euo pipefail

echo "Starting script execution..."

bash_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
echo "Script directory: ${bash_dir}"

# 复制文件
echo "Copying JE_model.m file..."
cp "${bash_dir}/JE_model.m" ./

if [ ! -f "JE_model.m" ]; then
    echo -e "${RED}Error: JE_model.m file not found after copying!${NC}"
    exit 1
fi

echo "Executing MATLAB script..."
# 运行MATLAB脚本
mrun JE_model.m
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to execute remake_onfly_den.m!${NC}"
    exit 1
fi

rm -f JE_model.m
echo "Temporary files removed."

echo -e "${GREEN}Successfully completed JE_model.m execution!${NC}"
echo "Script finished."