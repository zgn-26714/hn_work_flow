#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

int main(int argc, char *argv[]) {
    const char* v_str        = std::getenv("V");
    const char* time_str     = std::getenv("TIME");          // ps
    const char* dt_str       = std::getenv("DT");            // ps
    const char* skip_str     = std::getenv("SKIPTIME_PS");   // ps
    const char* tao_str      = std::getenv("TAO");           // "0 10 20"
    const char* outdir_str   = std::getenv("SLOWDIR");

    double V          = v_str    ? std::atof(v_str)    : 2.0;
    double total_time = time_str ? std::atof(time_str) : 100000.0; // ps
    double dt         = dt_str   ? std::atof(dt_str)   : 0.002;    // ps
    double skip_ps    = skip_str ? std::atof(skip_str) : 400.0;    // ps

    std::vector<double> tao;
    if (tao_str) {
        std::istringstream iss(tao_str);
        double v;
        while (iss >> v) tao.push_back(v);
    } else {
        tao = {0.0};
    }

    double Vini = 0.0;
    int skipTime = static_cast<int>(std::round(skip_ps / dt));
    int endT_val = static_cast<int>(std::round(total_time / dt));

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

        // 构建输出文件名
        std::string outfilename = std::string(outdir_str) + "Dphis_control.dat" + std::to_string(static_cast<int>(V)) 
                        + "_" + std::to_string(static_cast<int>(slow)) + "_" + std::to_string(std::stoi(skip_str));

        // 打开输出文件
        std::ofstream fout(outfilename);
        if (!fout.is_open()) {
            std::cerr << "Error opening file: " << outfilename << std::endl;
            continue;
        }

        // 写入文件头
        fout << "LINEAR\t# check if it's linear interpolate (more interpolate may be added?)\n";
        fout << (Vall.size()) + skipTime << "\t# total num of controlling points\n";

        // 设置浮点数输出格式（6位小数）
        fout << std::fixed << std::setprecision(6);
        long double timestamp = 0.0;
        fout<< timestamp<< " " << "0.000000" << "\n";
        // 写入数据点（跳过第一个时间点0）
        for(int i = 0; i < skipTime;i++){
            timestamp += dt; // 从0.002开始递增
            fout<< timestamp<< " " << "0.000000" << "\n";
        }
        

        for (int i = 1; i < Vall.size(); ++i) {
            timestamp += dt; // 从0.002开始递增
            fout << timestamp << " " << Vall[i] << "\n";
        }

        fout.close();
        std::cout << "done: " << outfilename << "\n";
    }

    return 0;
}
