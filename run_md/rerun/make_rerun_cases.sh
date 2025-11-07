#!/bin/bash

set -euo pipefail


# 假设这些函数在脚本的其他地方有定义
close_vdw() {
    local search_dir=$1
    local itp_files=()
    local default_files=()
    
    # Check if directory exists
    if [ ! -d "$search_dir" ]; then
        echo -e "${ERROR}Error: Directory '$search_dir' does not exist${NC}"
        return 1
    fi
    
    # Find all .itp files
    itp_files=($(find "$search_dir" -name "*.itp" -type f))
    
    if [ ${#itp_files[@]} -eq 0 ]; then
        echo -e "${ERROR}No .itp files found in directory '$search_dir'${NC}"
        return 1
    fi
    
    echo -e "${OK}Found ${#itp_files[@]} .itp files:"
    
    # Check each .itp file for "[ defaults ]" pattern
    for file in "${itp_files[@]}"; do
        if grep -q "\[ defaults \]" "$file" 2>/dev/null; then
            default_files+=("$file")
            echo -e "${OK}Found file with '[ defaults ]': $file"
        fi
    done
    
    if [ ${#default_files[@]} -eq 0 ]; then
        echo -e "${ERROR}No .itp files containing '[ defaults ]' found${NC}"
        return 1
    fi
    
    echo -e "${OK}Total files with '[ defaults ]' found: ${#default_files[@]}"
    
    # Store result in global variable
    CLOSE_VDW_RESULT=("${default_files[@]}")
    
    # 新增部分：如果只有一个文件，处理[ atomtypes ]部分
    if [ ${#CLOSE_VDW_RESULT[@]} -eq 1 ]; then
        local target_file="${CLOSE_VDW_RESULT[0]}"
        echo "Processing single file: $target_file"
        
        # 创建临时文件

        local temp_file=$(mktemp "${target_file}.XXXXXX")
        
        # 标志变量，用于跟踪是否在[ atomtypes ]部分
        local in_atomtypes_section=false
        
        # 读取文件并处理
        while IFS= read -r line; do
            # 检查是否进入[ atomtypes ]部分
            if [[ "$line" =~ ^\[[[:space:]]*atomtypes[[:space:]]*\] ]]; then
                in_atomtypes_section=true
                echo "$line"
                continue
            fi
            
            # 检查是否离开[ atomtypes ]部分（遇到下一个[ section ]）
            if [[ "$in_atomtypes_section" == true && "$line" =~ ^\[.*\] ]]; then
                in_atomtypes_section=false
            fi
            
            # 如果在[ atomtypes ]部分且不是注释行或空行
            if [[ "$in_atomtypes_section" == true && ! "$line" =~ ^[[:space:]]*\; && ! "$line" =~ ^[[:space:]]*$ ]]; then
                # 使用awk将最后两个数字清零
                local processed_line=$(echo "$line" | awk '{
                    # 计算字段数
                    nf = NF
                    # 如果最后两个字段是数字，将它们清零
                    if ($(nf-1) ~ /^[0-9.eE+-]+$/ && $nf ~ /^[0-9.eE+-]+$/) {
                        $(nf-1) = "0.00000"
                        $nf = "0.00000"
                    } else if ($(nf) ~ /^[0-9.eE+-]+;$/) {
                        # 处理以分号结尾的情况（如示例中的Cl行）
                        $(nf-1) = "0.00000"
                        $nf = "0.00000;"
                    }
                    # 打印处理后的行
                    for (i=1; i<=nf; i++) {
                        printf "%s", $i
                        if (i < nf) printf " "
                    }
                    printf "\n"
                }')
                echo "$processed_line"
            else
                echo "$line"
            fi
        done < "$target_file" > "$temp_file"
        
        # 用临时文件替换原文件
        mv "$temp_file" "$target_file"
        echo -e "${GREEN}Successfully processed [ atomtypes ] section in $target_file${NC}"
        
    elif [ ${#CLOSE_VDW_RESULT[@]} -gt 1 ]; then
        echo -e "${ERROR}Multiple files found with '[ defaults ]', skipping processing${NC}"
        echo "Files: ${CLOSE_VDW_RESULT[@]}"
    fi
}

close_bond() {
    local search_dir="$1"
    
    # 检查目录是否存在
    if [ ! -d "$search_dir" ]; then
        echo "Error: Directory '$search_dir' does not exist"
        return 1
    fi
    
    # 查找所有.itp文件
    local itp_files=($(find "$search_dir" -name "*.itp" -type f))
    
    if [ ${#itp_files[@]} -eq 0 ]; then
        echo "No .itp files found in directory '$search_dir'"
        return 1
    fi
    
    echo "Found ${#itp_files[@]} .itp files"
    
    # 处理每个文件
    for file in "${itp_files[@]}"; do
        echo "Processing: $file"
        
        # 先转换为Unix格式（去除\r）
        dos2unix "$file" 2>/dev/null || sed -i 's/\r$//' "$file" 2>/dev/null || true
        
        local temp_file="${file}.tmp"
        > "$temp_file"
        
        local in_section=""
        local keep_fields=0
        
        while IFS= read -r line || [ -n "$line" ]; do
            # 去除行尾的\r（如果dos2unix不可用）
            line="${line%$'\r'}"
            
            # 检测section
            if [[ "$line" =~ ^[[:space:]]*\[[[:space:]]*bonds[[:space:]]*\] ]]; then
                in_section="bonds"
                keep_fields=3  # ai, aj, funct
                echo "$line" >> "$temp_file"
                continue
            elif [[ "$line" =~ ^[[:space:]]*\[[[:space:]]*angles[[:space:]]*\] ]]; then
                in_section="angles"
                keep_fields=4  # ai, aj, ak, funct
                echo "$line" >> "$temp_file"
                continue
            elif [[ "$line" =~ ^[[:space:]]*\[[[:space:]]*dihedrals[[:space:]]*\] ]]; then
                in_section="dihedrals"
                keep_fields=5  # ai, aj, ak, al, funct
                echo "$line" >> "$temp_file"
                continue
            elif [[ "$line" =~ ^[[:space:]]*\[ ]]; then
                in_section=""
                keep_fields=0
                echo "$line" >> "$temp_file"
                continue
            fi
            
            # 如果不在目标section中，或者是注释/空行，直接输出
            if [[ -z "$in_section" ]] || [[ "$line" =~ ^[[:space:]]*\; ]] || [[ "$line" =~ ^[[:space:]]*$ ]]; then
                echo "$line" >> "$temp_file"
                continue
            fi
            
            # 处理数据行
            # 分离注释部分
            local data_part=""
            local comment_part=""
            
            if [[ "$line" =~ \; ]]; then
                data_part="${line%%;*}"
                comment_part=" ;${line##*;}"
            else
                data_part="$line"
                comment_part=""
            fi
            
            # 提取字段（去除多余空格）
            local fields=($data_part)
            local field_count=${#fields[@]}
            
            # 只有当字段数大于keep_fields时才需要修改
            if [ $field_count -le $keep_fields ]; then
                echo "$line" >> "$temp_file"
                continue
            fi
            
            # 构建输出：保留前keep_fields个字段，其余改为0
            local output=""
            for ((i=0; i<$field_count; i++)); do
                if [ $i -lt $keep_fields ]; then
                    output="$output ${fields[$i]}"
                else
                    output="$output 0"
                fi
            done
            
            # 输出处理后的行（去除开头空格后重新格式化）
            printf "%s%s\n" "$output" "$comment_part" >> "$temp_file"
            
        done < "$file"
        
        mv "$temp_file" "$file"
        echo "✓ Processed: $file"
    done
    
    echo "All files processed successfully!"
}

change_mdp() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "Error: file '$file' not found."
        return 1
    fi

    # 用 sed 逐行替换
    sed -i \
        -e 's/^[[:space:]]*nstxtcout[[:space:]]*=.*/nstxtcout       = 0/' \
        -e 's/^[[:space:]]*nstxout[[:space:]]*=.*/nstxout       = 1/' \
        -e 's/^[[:space:]]*nstfout[[:space:]]*=.*/nstfout       = 1/' \
        -e 's/^[[:space:]]*userint3[[:space:]]*=.*/userint3     = 1/' \
        -e 's/^[[:space:]]*userint2[[:space:]]*=.*/userint2    = 10000/' \
        "$file"
}



mkdir -p ./rerun_basicfile
if cp case${START}/*.top ./rerun_basicfile/topol.top; then
    echo "top file copied successfully"
else
    echo "top file copy failed, do you have more/less than one top in one case?"
    exit 1
fi
cp case${START}/*.itp ./rerun_basicfile
cp case${START}/*.ndx ./rerun_basicfile/index.ndx

if cp case${START}/grompp.mdp ./rerun_basicfile/case_rerun.mdp; then
    echo "mdp file copied successfully"
else
    echo "mdp file copy failed, please make sure you have grompp.mdp."
    exit 1
fi

close_vdw ./rerun_basicfile
close_bond ./rerun_basicfile
change_mdp ./rerun_basicfile/case_rerun.mdp

echo "当前路径: $(pwd)"

for (( i=$rerun_start; i<=$rerun_end; i++ )); do
    rm -rf ./case$i/rerun_case
    mkdir -p ./case$i/rerun_case
    cp ./rerun_basicfile/* ./case$i/rerun_case
    #CPM
    Matrix_dir=( $(realpath case$i/allMatrixA.bin) )
    control_file=( $(realpath case$i/CPM_ControlFile.dat) )
    Dphis=( $(realpath case$i/Dphis_control.dat) )
    ln -sf "${Matrix_dir}" case$i/rerun_case/allMatrixA.bin
    cp "${control_file}" case$i/rerun_case/CPM_ControlFile.dat
    ln -sf "${Dphis}" case$i/rerun_case/Dphis_control.dat
    ##traj
    shopt -s nullglob 

    traj_dir=( $(realpath case$i/*.xtc) )

    if [[ ${#traj_dir[@]} -eq 0 ]]; then
        echo -e "${YELLOW}⚠️  warning case$i! ln xtc failed, try to ln trr...${NC}"
        trr_traj_dir=(case$i/*.trr)
        if [[ ${#trr_traj_dir[@]} -eq 0 ]]; then
            echo "❌ Error: No .trr files found in case$i/"
            exit 1
        elif (( ${#trr_traj_dir[@]} > 1 )); then
            echo "❌ Error: Multiple .trr files found in case$i/:"
            printf '  %s\n' "${trr_traj_dir[@]}"
            exit 1
        else
            ln -sf "${trr_traj_dir[0]}" "./case$i/rerun_case/traj_to_rerun.trr"
            echo -e "${OK} trr for case$i" 
        fi
    elif (( ${#traj_dir[@]} > 1 )); then
        echo "❌ Error: Multiple .xtc files found in case$i/:"
        printf '  %s\n' "${traj_dir[@]}"
        exit 1
    else
        ln -sf "${traj_dir[0]}" "./case$i/rerun_case/traj_to_rerun.xtc"
        echo -e "${OK} xtc for case$i" 
    fi

done