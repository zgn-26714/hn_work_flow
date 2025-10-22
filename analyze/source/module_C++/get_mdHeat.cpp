#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "./include/matlab_to_C.h"

// 简单的二维数组结构来存储数据
struct Matrix {
    std::vector<std::vector<double>> data;
    
    size_t rows() const { return data.size(); }
    size_t cols() const { return data.empty() ? 0 : data[0].size(); }
    
    std::vector<double>& operator[](size_t index) { return data[index]; }
    const std::vector<double>& operator[](size_t index) const { return data[index]; }
};


int main() {
    const int begin_case = std::stoi(getenv("analyze_begin_case"));
    const int end_case = std::stoi(getenv("analyze_end_case"));
    
    std::vector<double> heat;
    Matrix data;
    bool first_case = true;
    
    int count = 0;
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
        data.data = matlab::xvgread(outedr);
        
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
                else {
                    std::cerr << "Warning: data is error! " << i << std::endl;
                    return 1;
                }
            }
            first_case = false;
        } else {
            for (size_t j = 0; j < data.rows() && j < heat.size(); ++j) {
                if (data[j].size() >= 3) {
                    heat[j] += data[j][1] - data[j][2];
                }
                else {
                    std::cerr << "Warning: data is error! " << i << std::endl;
                    return 1;//可以改掉
                }
            }
        }
        count++;
        // 清理临时文件
        system("rm #*# 2>/dev/null");
    }
    
    // 计算平均值
    for (size_t i = 0; i < heat.size(); ++i) {
        heat[i] /= count;
    }
    
    // 提取时间数据
    std::vector<double> time;
    if (data.rows() > 0) {
        for (size_t i = 0; i < data.rows(); ++i) {
            if (!data[i].empty()) {
                time.push_back(data[i][0]);  // 第一列是时间
            }
            else{
                std::cerr << "Warning: data is error when extracting time!" << std::endl;
                return 1;
            }
        }
    }
    

        // 保存为MAT文件
    std::string molName;
    std::string mat_filename = "./meandata/" + molName + "mdHeat" + getenv("V") + "V" + getenv("tau") + "ps.mat";
    if (matlab::saveMatFile(mat_filename, time, heat)) {
        std::cout << "MATLAB data saved to " <<mat_filename << std::endl;
    } else {
        std::cerr << "Failed to save MAT file" << std::endl;
    }
    
    return 0;
}