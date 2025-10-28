#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

// 原子结构体
struct Atom {
    int id;
    std::string name;
    std::string residue;
    int residue_id;
    double x, y, z;
    double occupancy;
    double temp_factor;
    double mass;  // 原子质量 (amu)
    
    // 获取原子电荷（简化版本，实际需要更精确的电荷参数）
    double getCharge() const {
        // 基于原子名称的简单电荷分配
        // 在实际应用中应该使用更精确的力场参数
        if (name == "C1" || name == "C2" || name == "C3" || name == "C4") {
            return 0.0;  // 碳原子，中性
        } else if (name[0] == 'O') {
            return -0.8; // 氧原子，带负电
        } else if (name[0] == 'H') {
            return 0.4;  // 氢原子，带正电
        }
        return 0.0; // 默认中性
    }
};

// 读取PDB文件
std::vector<Atom> readPDBFile(const std::string& filename) {
    std::vector<Atom> atoms;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return atoms;
    }
    
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM") {
            Atom atom;
            std::istringstream iss(line);
            std::string token;
            
            // 解析ATOM行
            iss >> token >> atom.id >> atom.name >> atom.residue >> atom.residue_id
                >> atom.x >> atom.y >> atom.z >> atom.occupancy >> atom.temp_factor;
            
            atom.mass = getAtomicMass(atom.name);
            atoms.push_back(atom);
        }
    }
    
    file.close();
    return atoms;
}

// 计算分子偶极矩
std::vector<double> calculateDipoleMoment(const std::vector<Atom>& atoms) {
    double dipole_x = 0.0, dipole_y = 0.0, dipole_z = 0.0;
    
    for (const auto& atom : atoms) {
        double charge = atom.getCharge();
        dipole_x += charge * atom.x /10;
        dipole_y += charge * atom.y /10;
        dipole_z += charge * atom.z /10;
    }
    
    return {dipole_x, dipole_y, dipole_z};
}

// 计算向量模长
double calculateMagnitude(const std::vector<double>& vector) {
    return std::sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}

// 打印结果
void printResults(const std::vector<double>& dipole, const std::vector<Atom>& atoms) {
    std::cout << "=== 分子偶极矩计算结果 ===" << std::endl;
    std::cout << "原子数量: " << atoms.size() << std::endl;
    std::cout << std::endl;
    
    std::cout << "各原子信息:" << std::endl;
    std::cout << "ID\tName\tRes\tX\t\tY\t\tZ\t\tCharge" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    
    for (const auto& atom : atoms) {
        printf("%d\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.2f\n", 
               atom.id, atom.name.c_str(), atom.residue.c_str(),
               atom.x, atom.y, atom.z, atom.getCharge());
    }
    
    std::cout << std::endl;
    std::cout << "偶极矩向量 (Debye):" << std::endl;
    printf("μ_x = %.4f D\n", dipole[0]);
    printf("μ_y = %.4f D\n", dipole[1]);
    printf("μ_z = %.4f D\n", dipole[2]);
    
    double magnitude = calculateMagnitude(dipole);
    printf("\n偶极矩大小: %.4f D\n", magnitude);
    
    // 计算单位方向向量
    if (magnitude > 1e-10) {
        std::vector<double> direction = {
            dipole[0] / magnitude,
            dipole[1] / magnitude,
            dipole[2] / magnitude
        };
        printf("方向向量: (%.4f, %.4f, %.4f)\n", direction[0], direction[1], direction[2]);
    }
}

int main() {
    std::vector<Atom> atoms = readPDBFile("xxx.pdb");
    std::vector<double> dipole = calculateDipoleMoment(atoms);

    std::vector<double> z_axis = {0, 0, 1};
    if(isparallel(dipole, z_axis)){
        std::cout << "parallel? or too small dipole" << std::endl;
        exit(1);
    } 

    std::vector<double> normal = crossProduct(dipole, z_axis);
    std::vector<double> unit_normal = normalize(normal);
    
    return 0;
}

std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2], 
        a[0]*b[1] - a[1]*b[0]
    };
}

bool isparallel(const std::vector<double>& a, const std::vector<double>& b) {
    double cross_x = a[1]*b[2] - a[2]*b[1];
    double cross_y = a[2]*b[0] - a[0]*b[2];
    double cross_z = a[0]*b[1] - a[1]*b[0];
    double magnitude = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
    return magnitude < 1e-10;
}

std::vector<double> normalize(const std::vector<double>& v) {
    double magnitude = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (magnitude < 1e-10) return {0, 0, 0};
    return {v[0]/magnitude, v[1]/magnitude, v[2]/magnitude};
}


double calculateMomentOfInertia(const std::vector<Atom>& atoms, 
                               const std::vector<double>& axis) {
    // 1. 计算质心
    std::vector<double> com = calculateCenterOfMass(atoms);
    
    // 2. 归一化旋转轴
    std::vector<double> unit_axis = normalize(axis);
    double ux = unit_axis[0], uy = unit_axis[1], uz = unit_axis[2];
    
    // 3. 计算转动惯量
    double I = 0.0;
    
    for (const auto& atom : atoms) {
        // 原子相对于质心的位置
        double rx = atom.x - com[0];
        double ry = atom.y - com[1];
        double rz = atom.z - com[2];
        
        // 到轴的距离的平方: r^2 - (r·u)^2
        double r_dot_u = rx*ux + ry*uy + rz*uz;
        double distance_sq = (rx*rx + ry*ry + rz*rz) - (r_dot_u * r_dot_u);
        
        I += atom.mass * distance_sq;
    }
    
    return I;
}


std::vector<double> calculateCenterOfMass(const std::vector<Atom>& atoms) {
    double total_mass = 0;
    double com_x = 0, com_y = 0, com_z = 0;
    
    for (const auto& atom : atoms) {
        com_x += atom.mass * atom.x;
        com_y += atom.mass * atom.y;
        com_z += atom.mass * atom.z;
        total_mass += atom.mass;
    }
    
    return {com_x/total_mass, com_y/total_mass, com_z/total_mass};
}

double getAtomicMass(const std::string& atom_name) {
    if (atom_name[0] == 'C') return 12.01;  // 碳
    if (atom_name[0] == 'O') return 16.00;  // 氧  
    if (atom_name[0] == 'H') return 1.008;  // 氢
    if (atom_name[0] == 'N') return 14.01;  // 氮
    if (atom_name[0] == 'P') return 30.97;  // 磷
    return 12.01; // 默认值
}