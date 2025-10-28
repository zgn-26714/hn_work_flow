#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <algorithm>

void merge_same_mol() {
    std::string top_file = getenv("TOP");
    std::string temp_file = top_file + ".tmp";
    
    std::cout << "开始合并top文件中相同分子的数量..." << std::endl;
    
    // 使用awk命令处理文件，合并相同分子
    std::string command = 
        "awk '"
        "BEGIN { in_molecules = 0 } "
        "/^\\[ molecules \\]/ { in_molecules = 1; print; next } "
        "/^\\[/ { in_molecules = 0 } "
        "in_molecules && /^[a-zA-Z][a-zA-Z0-9]*[ \\t]+[0-9]+/ { "
        "    molecule = $1; "
        "    count = $2; "
        "    if (molecule in counts) { "
        "        counts[molecule] += count "
        "    } else { "
        "        counts[molecule] = count "
        "    } "
        "    next "
        "} "
        "!in_molecules { print } "
        "END { "
        "    # 输出合并后的分子列表 "
        "    for (mol in counts) { "
        "        printf \"%-8s %d\\n\", mol, counts[mol] "
        "    } "
        "}' " + top_file + " > " + temp_file;
    
    // 执行合并操作
    int result = system(command.c_str());
    if (result == 0) {
        // 替换原文件
        system(("mv " + temp_file + " " + top_file).c_str());
        std::cout << "分子合并完成" << std::endl;
    } else {
        std::cerr << "合并操作失败" << std::endl;
        // 清理临时文件
        system(("rm -f " + temp_file).c_str());
    }
}



void change_top(const std::vector<std::string>& need_change_names, const std::vector<std::string>& target_names) {
    // 获取top文件名
    std::string top_file = getenv("TOP"); // 根据您的文件使用固定名称
    
    std::cout << "处理top文件: " << top_file << std::endl;
    
    // 检查文件是否存在
    std::ifstream file_check(top_file);
    if (!file_check.is_open()) {
        std::cerr << "错误：无法打开top文件 " << top_file << std::endl;
        return;
    }
    file_check.close();
    
    std::string temp_file = top_file + ".tmp";
    
    // 对每个需要替换的名称进行替换，使用与main函数相同的awk逻辑
    for (size_t i = 0; i < need_change_names.size(); ++i) {
        const std::string& old_name = need_change_names[i];
        const std::string& new_name = target_names[i];
        
        std::cout << "在top文件中替换: " << old_name << " -> " << new_name << std::endl;
        
        // 使用与main函数完全相同的awk命令模板
        std::string command = 
            "awk '{"
            "if (substr($0, 1, 20) ~ /[0-9]+" + old_name + "/) {"
            "pos = index($0, \"" + old_name + "\");"
            "if (pos > 0) {"
            "$0 = substr($0, 1, pos-1) \"" + new_name + "\" substr($0, pos+" + 
            std::to_string(old_name.length()) + ");"
            "}"
            "}"
            "print"
            "}' " + top_file + " > " + temp_file;
        
        // 执行替换
        system(command.c_str());
        system(("mv " + temp_file + " " + top_file).c_str());
    }
    
    std::cout << "top文件替换完成" << std::endl;
}


int main(int argc, char* argv[]) {
    // 定义需要替换的名称组
    std::vector<std::string> need_change_names;
    std::vector<std::string> target_names;

    std::string change_tmp = getenv("need_change_names");
    std::istringstream ss(change_tmp);
    std::string tmp;
    while(ss>>tmp){
        need_change_names.push_back(tmp);
    }
    
    std::string target_tmp = getenv("target_names");
    std::istringstream ss1(target_tmp);
    while(ss1>>tmp){
        target_names.push_back(tmp);
    }

    if (need_change_names.size() != target_names.size()) {
        std::cerr << "错误：两个向量长度不一致！" << std::endl;
        return 1;
    }
    
    std::string current_file = argv[1];
    std::string temp_file = current_file + ".tmp";
    
    for (size_t i = 0; i < need_change_names.size(); ++i) {
        const std::string& old_name = need_change_names[i];
        const std::string& new_name = target_names[i];
        
        std::cout << "替换: " << old_name << " -> " << new_name << std::endl;
        
        // 直接使用您的命令模板，只替换CL2和CL1部分
        std::string command = 
            "awk '{"
            "if (substr($0, 1, 20) ~ /[0-9]+" + old_name + "/) {"
            "pos = index($0, \"" + old_name + "\");"
            "if (pos > 0) {"
            "$0 = substr($0, 1, pos-1) \"" + new_name + "\" substr($0, pos+" + 
            std::to_string(old_name.length()) + ");"
            "}"
            "}"
            "print"
            "}' " + current_file + " > " + temp_file;
        
        system(command.c_str());
        system(("mv " + temp_file + " " + current_file).c_str());
    }
    
    change_top(need_change_names, target_names);
    merge_same_mol();
    std::cout << "all replace comlete." << std::endl;
    return 0;
}