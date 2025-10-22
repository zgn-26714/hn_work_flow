#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "matio.h"


namespace matlab{
    std::vector<std::vector<double>> xvgread(const std::string& filename, int lineex = 0) {
        std::vector<std::vector<double>> data;
        
        // Open file
        std::ifstream fid(filename);
        if (!fid.is_open()) {
            std::cout << "file not exist" << std::endl;
            return data;
        }
        
        int line = 0;
        int N = 0;
        std::string str;
        
        // First pass: find the starting line of data and number of columns
        while (std::getline(fid, str)) {
            line++;
            
            if (str.empty()) continue;
            
            char tmp = str[0];
            
            // Skip comment lines and special lines
            if (tmp != '#' && tmp != '@') {
                if (lineex == 0) {
                    // Count number of columns
                    std::istringstream iss(str);
                    double value;
                    N = 0;
                    while (iss >> value) {
                        N++;
                    }
                    break;
                } else {
                    lineex--;
                }
            }
        }
        
        fid.close();
        
        // Second pass: read data
        fid.open(filename);
        if (!fid.is_open()) {
            return data;
        }
        
        // Skip header lines
        int currentLine = 0;
        while (currentLine < line - 1 && std::getline(fid, str)) {
            currentLine++;
        }
        
        // Read data
        while (std::getline(fid, str)) {
            if (str.empty()) continue;
            
            char tmp = str[0];
            if (tmp == '#' || tmp == '@') continue;
            
            std::istringstream iss(str);
            std::vector<double> row;
            double value;
            
            while (iss >> value) {
                row.push_back(value);
            }
            
            if (!row.empty()) {
                data.push_back(row);
            }
        }
        
        fid.close();
        
        return data;
    }

    bool saveMatFile(const std::string& filename, 
                 const std::vector<double>& time, 
                 const std::vector<double>& heat) {
        mat_t *matfp;
        matfp = Mat_CreateVer(filename.c_str(), NULL, MAT_FT_MAT5);
        
        if (matfp == NULL) {
            std::cerr << "Error creating MAT file " << filename << std::endl;
            return false;
        }
        
        // 创建时间变量
        size_t dims_time[2] = {time.size(), 1};
        matvar_t *matvar_time = Mat_VarCreate("time", MAT_C_DOUBLE, MAT_T_DOUBLE, 
                                            2, dims_time, (void*)time.data(), 0);
        if (matvar_time != NULL) {
            Mat_VarWrite(matfp, matvar_time, MAT_COMPRESSION_ZLIB);
            Mat_VarFree(matvar_time);
        }
        
        // 创建heat变量
        size_t dims_heat[2] = {heat.size(), 1};
        matvar_t *matvar_heat = Mat_VarCreate("heat", MAT_C_DOUBLE, MAT_T_DOUBLE, 
                                            2, dims_heat, (void*)heat.data(), 0);
        if (matvar_heat != NULL) {
            Mat_VarWrite(matfp, matvar_heat, MAT_COMPRESSION_ZLIB);
            Mat_VarFree(matvar_heat);
        }
        
        Mat_Close(matfp);
        return true;
    }
}

