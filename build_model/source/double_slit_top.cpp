#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

int main() {
    std::string filename = getenv("TOP") + std::string(".top");
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法打开文件 '" << filename << "'" << std::endl;
        return 1;
    }

    std::vector<std::string> lines;
    std::string line;
    
    // 读取所有行
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();

    // 查找 [molecules] 段落开始位置
    int start_line = -1;
    for (size_t i = 0; i < lines.size(); ++i) {
        std::string trimmed = lines[i];
        trimmed.erase(std::remove(trimmed.begin(), trimmed.end(), ' '), trimmed.end());
        if (trimmed == "[molecules]") {
            start_line = i;
            break;
        }
    }

    if (start_line == -1) {
        std::cerr << "错误: 未找到 [molecules] 段落" << std::endl;
        return 1;
    }

    std::cout << "找到 molecules 段落: 从第 " << start_line + 1 << " 行到文件末尾" << std::endl;

    // 处理 molecules 段落
    std::vector<std::pair<std::string, int>> moleculeCounts;
    std::string firstMolecule;
    int firstCount = 0;
    bool firstMoleculeSet = false;

    // 从 molecules 段落开始处理到文件末尾
    for (size_t i = start_line + 1; i < lines.size(); ++i) {
        std::string currentLine = lines[i];
        
        // 跳过空行和注释行
        if (currentLine.empty() || currentLine[0] == ';') {
            continue;
        }

        // 跳过可能的新段落标题
        if (currentLine.find('[') != std::string::npos) {
            std::cerr << "error: find a [, check your file" << std::endl;
        }

        // 解析分子和数量
        std::istringstream iss(currentLine);
        std::string molecule;
        int count;
        
        if (iss >> molecule >> count) {
            std::pair<std::string, int> molPair = {molecule, count};
            auto iter = find(moleculeCounts.begin(), moleculeCounts.end(), molPair);
            if(iter == moleculeCounts.end()){
                moleculeCounts.push_back(molPair);
            }
            else{
                moleculeCounts[iter - moleculeCounts.begin()].second += count;
            }
            
            
            // 记录第一个分子
            if (!firstMoleculeSet) {
                firstMolecule = molecule;
                firstCount = count;
                firstMoleculeSet = true;
            }
        }
    }

    // 创建输出文件
    std::ofstream outFile(filename + ".processed");
    if (!outFile.is_open()) {
        std::cerr << "错误: 无法创建输出文件" << std::endl;
        return 1;
    }

    // 输出 molecules 段落之前的内容
    for (int i = 0; i < start_line; ++i) {
        outFile << lines[i] << std::endl;
    }

    // 输出 [molecules] 标题
    outFile << "[ molecules ]" << std::endl;

    // 输出处理后的 molecules 内容
    if (firstMoleculeSet) {
        // 输出第一个分子(EL)
        outFile << firstMolecule <<" "<< firstCount << std::endl;
        // 输出第一个分子的复制品(ER)
        outFile << "ER " << firstCount << std::endl;
    }

    // 输出其他分子，数量翻倍
    for (const auto& pair : moleculeCounts) {
        if (pair.first != firstMolecule) {
            outFile << pair.first << " " << pair.second * 2 << std::endl;
        }
    }

    outFile.close();
    std::cout << "处理完成！结果已保存到 " << filename << ".processed" << std::endl;

    return 0;
}