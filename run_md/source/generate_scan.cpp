#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstdlib> // for std::getenv, std::atof
#include <sstream> // for std::istringstream
#include <filesystem>

namespace fs = std::filesystem;
int main() {
    // 从环境变量获取值，如果不存在则使用默认值
    const char* v_str = std::getenv("V");
    const char* time_str = std::getenv("TIME");
    const char* dt_str = std::getenv("DT");
    const char* skipTime_str = std::getenv("SKIPTIME_PS");
    const char* tao_str = std::getenv("TAO");
    const char* slowdir = std::getenv("SLOWDIR");

    double V = v_str ? std::atof(v_str) : 2.0;
    double time = time_str ? std::atof(time_str) : 10000.0; // ps
    double dt = dt_str ? std::atof(dt_str) : 0.002; // ps
    double skipTime_ps = skipTime_str ? std::atof(skipTime_str) : 400.0; // ps
    
    // 处理tao值（空格分隔的多个值）
    std::vector<double> tao;
    if (tao_str) {
        std::string tao_env = tao_str;
        std::istringstream iss(tao_env);
        double value;
        while (iss >> value) {
            tao.push_back(value);
        }
    } else {
        tao = {0}; // 默认值
    }

    double Vini = 0.0;
    double endT_val = time / dt;
    int skipTime = static_cast<int>(std::round(skipTime_ps / dt));
    

    for (double slow : tao) {
        double tini = 0.0;
        double tup = slow;
        double tend = 0.0;

        double tall = tup + tini + tend;
        
        // 计算时间点数量（整数）
        int n = static_cast<int>(std::round(tall / dt)) + 1;
        std::vector<double> Vall(n, 0.0);

        // 生成电压序列
        for (int i = 0; i < n; ++i) {
            double nowt = i * dt;
            if (nowt < tini) {
                Vall[i] = Vini;
            } else if (nowt > tini + tup) {
                Vall[i] = V;
            } else {
                Vall[i] = Vini + (nowt - tini) / tup * (V - Vini);
            }
        }

        // 检查并补全电压序列
        if (Vall.size() < endT_val) {
            Vall.resize(static_cast<int>(endT_val), V);
        }

        std::string outfilename = std::string(slowdir) + "Dphis_control.dat" + std::to_string(static_cast<int>(V)) 
                    + "_" + std::to_string(static_cast<int>(slow)) + "_" + std::to_string(std::stoi(skipTime_str));
        // 打开输出文件
        if (fs::exists(outfilename)) {
            std::cout << "file exist " << outfilename << ", continue..." << std::endl;
            continue;
        }
        std::cout << "file not exist " << outfilename << ", generate..." << std::endl;

        std::ofstream fout(outfilename);
        if (!fout.is_open()) {
            std::cerr << "Error opening file: " << outfilename << std::endl;
            continue;
        }
	    std::cout<<"max simulation step:\t"<<(Vall.size() - 1)<<'\n'<<"skipstep:\t"<<skipTime<<std::endl;
        // 写入文件头
        fout << "numofdphis\t";
        fout << (Vall.size() - 1) + skipTime<<std::endl;

        // 设置浮点数输出格式（6位小数）
        fout << std::fixed << std::setprecision(6);
        // 写入数据点（跳过第一个时间点0）
        for(int i = 0; i < skipTime;i++){
            fout<< "0.000000" << "\n";
        }
        

        for (int i = 1; i < Vall.size(); ++i) {
            fout << Vall[i] << "\n";
        }

        fout.close();
    }

    return 0;
}
