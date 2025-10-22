#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>

// 简单的二维数组结构来存储数据
struct Matrix {
    std::vector<std::vector<double>> data;
    
    size_t rows() const { return data.size(); }
    size_t cols() const { return data.empty() ? 0 : data[0].size(); }
    
    std::vector<double>& operator[](size_t index) { return data[index]; }
    const std::vector<double>& operator[](size_t index) const { return data[index]; }
};

// 模拟MATLAB的xvgread函数
Matrix xvgread(const std::string& filename) {
    Matrix result;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return result;
    }
    
    // 跳过注释行（以#、@开头的行）
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '@') {
            continue;
        }
        
        std::vector<double> row;
        std::istringstream iss(line);
        double value;
        
        while (iss >> value) {
            row.push_back(value);
        }
        
        if (!row.empty()) {
            result.data.push_back(row);
        }
    }
    
    file.close();
    return result;
}

int main() {
    const int begin_case = std::stoi(getenv("analyze_begin_case"));
    const int end_case = std::stoi(getenv("analyze_end_case"));
    
    std::vector<double> heat;
    Matrix data;
    bool first_case = true;
    
    for (int i = begin_case; i <= end_case; ++i) {
        // 构建edr文件路径
        
        std::ostringstream edr_oss;
        std::string deffnm = getenv("DEFFNM");
        edr_oss << "./case" << i << "/" << deffnm << ".edr";
        std::string edrdir = edr_oss.str();
        
        std::string outedr = "./energy.xvg";
        
        // 构建并执行GROMACS命令
        std::ostringstream command_oss;
        command_oss << "echo -e \"Total-Energy\\nConserved-En\" | gmx energy -f " 
                   << edrdir << " -o " << outedr;
        
        int result = system(command_oss.str().c_str());
        
        if (result != 0) {
            std::cerr << "Warning: GROMACS command failed for case " << i << std::endl;
            continue;
        }
        
        // 读取xvg文件
        data = xvgread(outedr);
        
        if (data.rows() == 0) {
            std::cerr << "Warning: No data read for case " << i << std::endl;
            continue;
        }
        
        // 处理数据
        if (first_case) {
            heat.resize(data.rows());
            for (size_t j = 0; j < data.rows(); ++j) {
                if (data[j].size() >= 3) {
                    heat[j] = data[j][1] - data[j][2];  // 第2列 - 第3列
                }
            }
            first_case = false;
        } else {
            for (size_t j = 0; j < data.rows() && j < heat.size(); ++j) {
                if (data[j].size() >= 3) {
                    heat[j] += data[j][1] - data[j][2];
                }
            }
        }
        
        // 清理临时文件
        system("rm *~ *# 2>/dev/null");
    }
    
    // 计算平均值
    int num_cases = end_case - begin_case + 1;
    for (size_t i = 0; i < heat.size(); ++i) {
        heat[i] /= num_cases;
    }
    
    // 提取时间数据
    std::vector<double> time;
    if (data.rows() > 0) {
        for (size_t i = 0; i < data.rows(); ++i) {
            if (!data[i].empty()) {
                time.push_back(data[i][0]);  // 第一列是时间
            }
        }
    }
    
    // 保存数据到文件（简化版本）
    std::ofstream outfile("./meandata/ACNHeat2V0ps.txt");
    if (outfile.is_open()) {
        outfile << "Time\tHeat" << std::endl;
        for (size_t i = 0; i < time.size() && i < heat.size(); ++i) {
            outfile << std::scientific << std::setprecision(6) 
                   << time[i] << "\t" << heat[i] << std::endl;
        }
        outfile.close();
    }
    
    // 如果需要MATLAB格式，可以添加MAT文件输出
    // 这里需要额外的库如matio来写入.mat文件
    
    std::cout << "Processing completed. Results saved to ./meandata/ACNHeat2V0ps.txt" << std::endl;
    
    return 0;
}