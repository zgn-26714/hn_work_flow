#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>

// 原子坐标结构体
struct Atom {
    int resNum;           // 残基编号
    std::string resName;  // 残基名称
    std::string atomName; // 原子名称
    int atomNum;          // 原子编号
    double x, y, z;       // 坐标 (nm)
};

// 模拟盒子
struct Box {
    double x, y, z;
};

// Nlist结构：存储不同类型原子的数量
struct NList {
    std::vector<int> counts;
};

// 从文件读取单元结构
std::vector<std::vector<double>> loadUnitFile(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "警告: 无法打开文件 " << filename << std::endl;
        return data;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double val;
        while (iss >> val) {
            row.push_back(val);
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    file.close();
    return data;
}

// 生成石墨烯电极结构
void GenElec_zl(double llx, double llz, std::vector<Atom>& Coor, Box& box, NList& Nlist) {
    // 加载单元结构文件
    auto Unit1 = loadUnitFile("./database/Unit1_cvo.dat");
    auto Unit2 = loadUnitFile("./database/Unit2_cv.dat");
    auto Unit3 = loadUnitFile("./database/Unit3_cv.dat");
    auto Unit4 = loadUnitFile("./database/Unit4_cpo.dat");
    auto Unit5 = loadUnitFile("./database/Unit5_cp.dat");
    auto Unit6 = loadUnitFile("./database/Unit6_cp.dat");
    
    // 参数设置
    const double Lx0 = 3.3;
    const double Ly0 = 3.0;
    const double Lz0 = 6.1;
    const double dc = 0.2460;
    const double dg = 0.3410;
    
    // 计算网格扩展
    int nx = static_cast<int>(std::round((llx - Lx0) / dc + 7));
    int nz = static_cast<int>(std::round((llz - Lz0) / dc + 19));
    
    // 计算各类型原子数量
    int N1 = 14 * (2 * (13 + nx - 7) + 1);
    int N2 = 14 * (2 * nx + 1);
    int N3 = 14 * (2 * (10 + nx - 7) + 1);
    int N4 = 14 * (2 * (24 + nz - 19) + 1);
    int N5 = 14 * (2 * (21 + nz - 19) + 1);
    int N6 = 14 * (2 * nz + 1);
    
    int Natom = 2 * (N1 + N2 + N3 + N4 + N5 + N6);
    
    // 设置Nlist
    Nlist.counts = {2*N4, 2*N1, 2*(N5+N6), 2*(N3+N2)};
    
    // 设置盒子大小
    box.x = Lx0 + (nx - 7) * dc;
    box.y = Ly0;
    box.z = Lz0 + (nz - 19) * dc;
    
    Coor.clear();
    Coor.resize(Natom);
    
    int idx = 0;
    
    // ===== Cpo 区域 =====
    // MATLAB: for i=1:24+nz-19, for j=1:28
    // C++: i从0开始，对应MATLAB的i=1; j从0开始对应j=1
    // 第一部分
    for (int i = 1; i <= 24 + nz - 19; i++) {
        for (int j = 1; j <= 28; j++) {
            if (j-1 < Unit4.size() && idx < Natom) {
                Coor[idx].resNum = idx + 1;
                Coor[idx].resName = "Cpo";
                Coor[idx].atomName = "C";
                Coor[idx].atomNum = idx + 1;
                Coor[idx].x = Unit4[j-1][0];  // MATLAB j对应C++ j-1
                Coor[idx].y = Unit4[j-1][1];
                Coor[idx].z = Unit4[j-1][2] + (i-1) * dc;  // MATLAB i对应C++ i-1
                idx++;
            }
        }
    }
    
    // 第二部分：MATLAB for i=1:14
    for (int i = 1; i <= 14 && i-1 < Unit4.size(); i++) {
        if (idx < Natom) {
            Coor[idx].resNum = idx + 1;
            Coor[idx].resName = "Cpo";
            Coor[idx].atomName = "C";
            Coor[idx].atomNum = idx + 1;
            Coor[idx].x = Unit4[i-1][0];  // MATLAB i对应C++ i-1
            Coor[idx].y = Unit4[i-1][1];
            Coor[idx].z = Unit4[i-1][2] + (24 + nz - 19) * dc;
            idx++;
        }
    }
    
    // 镜像部分：MATLAB for i=1:14*(2*(24+nz-19)+1)
    int startCpo = 0;
    for (int i = 1; i <= N4 && idx < Natom; i++) {
        Coor[idx].resNum = idx + 1;
        Coor[idx].resName = "Cpo";
        Coor[idx].atomName = "C";
        Coor[idx].atomNum = idx + 1;
        Coor[idx].x = Lx0 + (nx - 7) * dc;
        Coor[idx].y = Coor[startCpo + i - 1].y;  // MATLAB的i对应C++的i-1
        Coor[idx].z = Coor[startCpo + i - 1].z;
        idx++;
    }
    
    // ===== Cvo 区域 =====
    int startCvo = idx;
    // MATLAB: for i=1:13+nx-7, for j=1:28
    for (int i = 1; i <= 13 + nx - 7; i++) {
        for (int j = 1; j <= 28; j++) {
            if (j-1 < Unit1.size() && idx < Natom) {
                Coor[idx].resNum = idx + 1;
                Coor[idx].resName = "Cvo";
                Coor[idx].atomName = "C";
                Coor[idx].atomNum = idx + 1;
                Coor[idx].x = Unit1[j-1][0] + (i-1) * dc;
                Coor[idx].y = Unit1[j-1][1];
                Coor[idx].z = Unit1[j-1][2];
                idx++;
            }
        }
    }
    
    // MATLAB: for i=1:14
    for (int i = 1; i <= 14 && i-1 < Unit1.size(); i++) {
        if (idx < Natom) {
            Coor[idx].resNum = idx + 1;
            Coor[idx].resName = "Cvo";
            Coor[idx].atomName = "C";
            Coor[idx].atomNum = idx + 1;
            Coor[idx].x = Unit1[i-1][0] + (13 + nx - 7) * dc;
            Coor[idx].y = Unit1[i-1][1];
            Coor[idx].z = Unit1[i-1][2];
            idx++;
        }
    }
    
    // MATLAB: for i=1:14*(2*(13+nx-7)+1)
    for (int i = 1; i <= N1 && idx < Natom; i++) {
        Coor[idx].resNum = idx + 1;
        Coor[idx].resName = "Cvo";
        Coor[idx].atomName = "C";
        Coor[idx].atomNum = idx + 1;
        Coor[idx].x = Coor[startCvo + i - 1].x;
        Coor[idx].y = Coor[startCvo + i - 1].y;
        Coor[idx].z = Lz0 + (nz - 19) * dc;
        idx++;
    }
    
    // ===== Cp 区域 (Unit5) =====
    int startCp5 = idx;
    // MATLAB: for i=1:21+nz-19
    for (int i = 1; i <= 21 + nz - 19; i++) {
        for (int j = 1; j <= 28; j++) {
            if (j-1 < Unit5.size() && idx < Natom) {
                Coor[idx].resNum = idx + 1;
                Coor[idx].resName = "Cp";
                Coor[idx].atomName = "C";
                Coor[idx].atomNum = idx + 1;
                Coor[idx].x = Unit5[j-1][0];
                Coor[idx].y = Unit5[j-1][1];
                Coor[idx].z = Unit5[j-1][2] + (i-1) * dc;
                idx++;
            }
        }
    }
    
    // MATLAB: for i=1:14
    for (int i = 1; i <= 14 && i-1 < Unit5.size(); i++) {
        if (idx < Natom) {
            Coor[idx].resNum = idx + 1;
            Coor[idx].resName = "Cp";
            Coor[idx].atomName = "C";
            Coor[idx].atomNum = idx + 1;
            Coor[idx].x = Unit5[i-1][0];
            Coor[idx].y = Unit5[i-1][1];
            Coor[idx].z = Unit5[i-1][2] + (21 + nz - 19) * dc;
            idx++;
        }
    }
    
    // MATLAB: for i=1:14*(2*(21+nz-19)+1)
    for (int i = 1; i <= N5 && idx < Natom; i++) {
        Coor[idx].resNum = idx + 1;
        Coor[idx].resName = "Cp";
        Coor[idx].atomName = "C";
        Coor[idx].atomNum = idx + 1;
        Coor[idx].x = Lx0 + (nx - 7) * dc - dg;
        Coor[idx].y = Coor[startCp5 + i - 1].y;
        Coor[idx].z = Coor[startCp5 + i - 1].z;
        idx++;
    }
    
    // ===== Cp 区域 (Unit6) =====
    int startCp6 = idx;
    // MATLAB: for i=1:19+nz-19
    for (int i = 1; i <= 19 + nz - 19; i++) {
        for (int j = 1; j <= 28; j++) {
            if (j-1 < Unit6.size() && idx < Natom) {
                Coor[idx].resNum = idx + 1;
                Coor[idx].resName = "Cp";
                Coor[idx].atomName = "C";
                Coor[idx].atomNum = idx + 1;
                Coor[idx].x = Unit6[j-1][0];
                Coor[idx].y = Unit6[j-1][1];
                Coor[idx].z = Unit6[j-1][2] + (i-1) * dc;
                idx++;
            }
        }
    }
    
    // MATLAB: for i=1:14
    for (int i = 1; i <= 14 && i-1 < Unit6.size(); i++) {
        if (idx < Natom) {
            Coor[idx].resNum = idx + 1;
            Coor[idx].resName = "Cp";
            Coor[idx].atomName = "C";
            Coor[idx].atomNum = idx + 1;
            Coor[idx].x = Unit6[i-1][0];
            Coor[idx].y = Unit6[i-1][1];
            Coor[idx].z = Unit6[i-1][2] + nz * dc;
            idx++;
        }
    }
    
    // MATLAB: for i=1:14*(2*nz+1)
    for (int i = 1; i <= N6 && idx < Natom; i++) {
        Coor[idx].resNum = idx + 1;
        Coor[idx].resName = "Cp";
        Coor[idx].atomName = "C";
        Coor[idx].atomNum = idx + 1;
        Coor[idx].x = Lx0 + (nx - 7) * dc - 2 * dg;
        Coor[idx].y = Coor[startCp6 + i - 1].y;
        Coor[idx].z = Coor[startCp6 + i - 1].z;
        idx++;
    }
    
    // ===== Cv 区域 (Unit3) =====
    int startCv3 = idx;
    // MATLAB: for i=1:10+nx-7
    for (int i = 1; i <= 10 + nx - 7; i++) {
        for (int j = 1; j <= 28; j++) {
            if (j-1 < Unit3.size() && idx < Natom) {
                Coor[idx].resNum = idx + 1;
                Coor[idx].resName = "Cv";
                Coor[idx].atomName = "C";
                Coor[idx].atomNum = idx + 1;
                Coor[idx].x = Unit3[j-1][0] + (i-1) * dc;
                Coor[idx].y = Unit3[j-1][1];
                Coor[idx].z = Unit3[j-1][2];
                idx++;
            }
        }
    }
    
    // MATLAB: for i=1:14
    for (int i = 1; i <= 14 && i-1 < Unit3.size(); i++) {
        if (idx < Natom) {
            Coor[idx].resNum = idx + 1;
            Coor[idx].resName = "Cv";
            Coor[idx].atomName = "C";
            Coor[idx].atomNum = idx + 1;
            Coor[idx].x = Unit3[i-1][0] + (10 + nx - 7) * dc;
            Coor[idx].y = Unit3[i-1][1];
            Coor[idx].z = Unit3[i-1][2];
            idx++;
        }
    }
    
    // MATLAB: for i=1:14*(2*(10+nx-7)+1)
    for (int i = 1; i <= N3 && idx < Natom; i++) {
        Coor[idx].resNum = idx + 1;
        Coor[idx].resName = "Cv";
        Coor[idx].atomName = "C";
        Coor[idx].atomNum = idx + 1;
        Coor[idx].x = Coor[startCv3 + i - 1].x;
        Coor[idx].y = Coor[startCv3 + i - 1].y;
        Coor[idx].z = Lz0 + (nz - 19) * dc - dg;
        idx++;
    }
    
    // ===== Cv 区域 (Unit2) =====
    int startCv2 = idx;
    // MATLAB: for i=1:7+nx-7
    for (int i = 1; i <= 7 + nx - 7; i++) {
        for (int j = 1; j <= 28; j++) {
            if (j-1 < Unit2.size() && idx < Natom) {
                Coor[idx].resNum = idx + 1;
                Coor[idx].resName = "Cv";
                Coor[idx].atomName = "C";
                Coor[idx].atomNum = idx + 1;
                Coor[idx].x = Unit2[j-1][0] + (i-1) * dc;
                Coor[idx].y = Unit2[j-1][1];
                Coor[idx].z = Unit2[j-1][2];
                idx++;
            }
        }
    }
    
    // MATLAB: for i=1:14
    for (int i = 1; i <= 14 && i-1 < Unit2.size(); i++) {
        if (idx < Natom) {
            Coor[idx].resNum = idx + 1;
            Coor[idx].resName = "Cv";
            Coor[idx].atomName = "C";
            Coor[idx].atomNum = idx + 1;
            Coor[idx].x = Unit2[i-1][0] + nx * dc;
            Coor[idx].y = Unit2[i-1][1];
            Coor[idx].z = Unit2[i-1][2];
            idx++;
        }
    }
    
    // MATLAB: for i=1:14*(2*nx+1)
    for (int i = 1; i <= N2 && idx < Natom; i++) {
        Coor[idx].resNum = idx + 1;
        Coor[idx].resName = "Cv";
        Coor[idx].atomName = "C";
        Coor[idx].atomNum = idx + 1;
        Coor[idx].x = Coor[startCv2 + i - 1].x;
        Coor[idx].y = Coor[startCv2 + i - 1].y;
        Coor[idx].z = Lz0 + (nz - 19) * dc - 2 * dg;
        idx++;
    }
    
    // 调整实际生成的原子数
    Coor.resize(idx);
}

// 平移并处理周期性边界条件
void translate_pbc(const std::vector<Atom>& Coor_sq, const Box& box_sq, double dslit,
                   std::vector<Atom>& Coor, Box& box) {
    double bulk_z = box_sq.z;
    Coor = Coor_sq;
    
    // 平移z坐标
    for (auto& atom : Coor) {
        atom.z = atom.z + bulk_z / 2.0;
    }
    
    // 平移x坐标
    for (auto& atom : Coor) {
        atom.x = atom.x + box_sq.x / 2.0 + dslit;
    }
    
    // 设置新盒子尺寸
    box.x = box_sq.x + dslit;
    box.y = box_sq.y;
    box.z = 2 * box_sq.z;
    
    // 处理周期性边界条件（仅对x方向）
    for (auto& atom : Coor) {
        if (atom.x > box.x) {
            atom.x = atom.x - box.x;
        }
    }
}

// 简化命名
void simplename(const std::vector<Atom>& Coor2, const Box& box2, const NList& Nlist,
                std::vector<Atom>& Coor4, Box& box4) {
    Coor4 = Coor2;
    box4 = box2;
    
    std::vector<std::string> Mname = {"CL2", "CR2", "CL1", "CR1"};
    
    int N = Coor2.size();
    int N1 = Nlist.counts[0];
    int sumN = 0;
    for (int count : Nlist.counts) {
        sumN += count;
    }
    
    // 重命名不同区域的原子
    // 区域1: 0 到 N1-1
    for (int i = 0; i < N1 && i < Coor4.size(); i++) {
        Coor4[i].resName = Mname[0];
    }
    
    // 区域2: N1 到 2*N1-1
    for (int i = N1; i < 2 * N1 && i < Coor4.size(); i++) {
        Coor4[i].resName = Mname[1];
    }
    
    // 区域3: 2*N1 到 N1+sumN-1
    for (int i = 2 * N1; i < N1 + sumN && i < Coor4.size(); i++) {
        Coor4[i].resName = Mname[2];
    }
    
    // 区域4: N1+sumN 到 2*sumN-1
    for (int i = N1 + sumN; i < 2 * sumN && i < Coor4.size(); i++) {
        Coor4[i].resName = Mname[3];
    }
}

// 写入.gro文件
void growrite(const std::vector<Atom>& Coor, const Box& box, const std::string& fname) {
    std::ofstream fout(fname);
    if (!fout.is_open()) {
        std::cerr << "无法打开文件: " << fname << std::endl;
        return;
    }
    
    fout << "GrapheneElectrodeSingle" << std::endl;
    fout << std::setw(5) << Coor.size() << std::endl;
    
    for (const auto& atom : Coor) {
        fout << std::setw(5) << atom.resNum
             << std::left << std::setw(5) << atom.resName
             << std::right << std::setw(5) << atom.atomName
             << std::setw(5) << atom.atomNum
             << std::fixed << std::setprecision(3)
             << std::setw(8) << atom.x
             << std::setw(8) << atom.y
             << std::setw(8) << atom.z
             << std::endl;
    }
    
    fout << std::fixed << std::setprecision(5)
         << std::setw(10) << box.x
         << std::setw(10) << box.y
         << std::setw(10) << box.z
         << std::endl;
    
    fout.close();
}

// 主函数示例
int main() {
    double llx = 3.4;  // nm
    double llz = 8.0;  // nm
    double dslit = 0.8; // nm
    
    std::vector<Atom> Coor_sq, Coor;
    Box box_sq, box;
    NList Nlist;
    
    // 生成初始电极结构
    GenElec_zl(llx, llz, Coor_sq, box_sq, Nlist);
    
    // 添加狭缝并平移
    translate_pbc(Coor_sq, box_sq, dslit, Coor, box);
    
    // 写入文件
    growrite(Coor, box, "oneelectrode.gro");
    
    std::cout << "生成了 " << Coor.size() << " 个原子" << std::endl;
    std::cout << "盒子尺寸: " << box.x << " x " << box.y << " x " << box.z << " nm^3" << std::endl;
    
    return 0;
}