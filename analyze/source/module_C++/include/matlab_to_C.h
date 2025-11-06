#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>


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
        
        // // First pass: find the starting line of data and number of columns
        // while (std::getline(fid, str)) {
        //     line++;
            
        //     if (str.empty()) continue;
            
        //     char tmp = str[0];
            
        //     // Skip comment lines and special lines
        //     if (tmp != '#' && tmp != '@') {
        //         if (lineex == 0) {
        //             // Count number of columns
        //             std::istringstream iss(str);
        //             double value;
        //             N = 0;
        //             while (iss >> value) {
        //                 N++;
        //             }
        //             break;
        //         } else {
        //             lineex--;
        //         }
        //     }
        // }
        
        // fid.close();
        
        // // Second pass: read data
        // fid.open(filename);
        // if (!fid.is_open()) {
        //     return data;
        // }
        
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

    bool xvgsave(const std::string& filename, 
                 const std::vector<std::string>& header, 
                 const std::vector<std::vector<double>>& data) {
        
        std::ofstream fid(filename);
        if (!fid.is_open()) {
            return false;
        }
        
        // Write header
        for (const auto& line : header) {
            fid << line << std::endl;
        }
        
        // Write data
        for (const auto& row : data) {
            for (size_t j = 0; j < row.size(); ++j) {
                fid << row[j];
                if (j < row.size() - 1) {
                    fid << " ";
                }
            }
            fid << std::endl;
        }
        
        fid.close();
        return true;
    }


    bool saveDatFile(const std::string& filename, 
                 const std::vector<double>& time, 
                 const std::vector<double>& heat) {
        
        return true;
    }
}

