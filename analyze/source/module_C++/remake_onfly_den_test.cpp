#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <regex>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

// 4D数组: [nx][ny][nz][col]
using Matrix4D = vector<vector<vector<vector<double>>>>;
using Matrix3D = vector<vector<vector<double>>>;

class remake_data
{
struct mol_data
{
    vector<vector<double>> zt;
    Matrix3D xyz;
    Matrix4D allData;

    vector<vector<double>>& get_Matrix(int dim1, int dim2) {
        zt.resize(dim1, vector<double>(dim2, 0.0));
        return zt;
    }

    Matrix3D& get_xyz(int dim1, int dim2, int dim3) {
        xyz.resize(dim1, vector<vector<double>>(dim2, vector<double>(dim3, 0.0)));
        return xyz;
    }

    Matrix4D& get_allData(int dim1, int dim2, int dim3, int dim4) {
        allData.resize(dim1, vector<vector<vector<double>>>(dim2, vector<vector<double>>(dim3, vector<double>(dim4, 0.0))));
        return allData;
    }

    mol_data& operator=(const mol_data& other) {
        if (this != &other) {
            zt = other.zt;
            xyz = other.xyz;
            allData = other.allData;
        }
        return *this;
    }

    mol_data& operator+=(const mol_data& other){
        if (this != &other) {
            // Assuming same dimensions for simplicity
            for (size_t i = 0; i < zt.size(); ++i) {
                for (size_t j = 0; j < zt[i].size(); ++j) {
                    zt[i][j] += other.zt[i][j];
                }
            }
            // Similar addition can be implemented for xyz and allData if needed
            for (size_t i = 0; i < xyz.size(); ++i) {
                for (size_t j = 0; j < xyz[i].size(); ++j) {
                    for (size_t k = 0; k < xyz[i][j].size(); ++k) {
                        xyz[i][j][k] += other.xyz[i][j][k];
                    }
                }
            }
            for (size_t i = 0; i < allData.size(); ++i) {
                for (size_t j = 0; j < allData[i].size(); ++j) {
                    for (size_t k = 0; k < allData[i][j].size(); ++k) {
                        for (size_t l = 0; l < allData[i][j][k].size(); ++l) {
                            allData[i][j][k][l] += other.allData[i][j][k][l];
                        }
                    }
                }
            }
        }
        return *this;
    }
};

public:
    vector<vector<mol_data>> data;
    remake_data(int nRegion, int nMolecule) {
        data.resize(nRegion, vector<mol_data>(nMolecule));
    }
    vector<mol_data>& get_region(int region) {
        return data[region]; 
    }
};



struct LoadDataResult {
    vector<vector<Matrix4D>> ret;
    vector<vector<double>> xyzRange;
    vector<int> nbin;
    vector<double> lowPos;
    vector<double> upPos;
    vector<double> Lbox;
    int totNbin;
    int totFrames;
    int nRegion;
    int nGroups;
};

vector<double> extractNumbers(const string& str) {
    vector<double> numbers;
    regex number_regex(R"([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)");
    sregex_iterator iter(str.begin(), str.end(), number_regex);
    sregex_iterator end;
    for (; iter != end; ++iter) {
        numbers.push_back(stod(iter->str()));
    }
    return numbers;
}

vector<string> readTextData(const string& filename) {
    vector<string> textData;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        if (!line.empty() && (line[0] == '%' || line[0] == '#')) {
            textData.push_back(line);
        }
        else {
            break;  // Stop reading on first non-comment line,it's mean you shoudn't modified the data file format
        }
    }
    file.close();
    return textData;
}

vector<vector<double>> readDataFile(const string& filename) {
    vector<vector<double>> data;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '%' || line[0] == '#') continue;
        
        vector<double> row;
        istringstream iss(line);
        double value;
        while (iss >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    file.close();
    return data;
}

vector<Matrix4D> reshapeColumnMajor(
    const vector<vector<double>>& tmpData,
    int nRegion,
    const vector<int>& nbin,
    int col)
{
    vector<Matrix4D> ret;
    ret.reserve(nRegion);

    size_t tmpIndex = 0;
    size_t totalRows = tmpData.size();

    for (int ii = 0; ii < nRegion; ++ii)
    {
        int n1 = nbin[3 * ii];       // 对应 nbinX
        int n2 = nbin[3 * ii + 1];   // 对应 nbinY
        int n3 = nbin[3 * ii + 2];   // 对应 nbinZ
        size_t eachNbin = static_cast<size_t>(n1) * n2 * n3;

        if (tmpIndex + eachNbin > totalRows)
            throw runtime_error("tmpData size not sufficient for nbin region");

        // 分配 4D 容器
        Matrix4D regionData(
            n1, vector<vector<vector<double>>>(
            n2, vector<vector<double>>(
            n3, vector<double>(col, 0.0))));

        // MATLAB 是列主序: x 最快，y 次，z 次，col 最慢
        size_t idx = 0;
        for (int c = 0; c < col; ++c){
            for (int z = 0; z < n3; ++z)
                for (int y = 0; y < n2; ++y)
                    for (int x = 0; x < n1; ++x)
                    {
                        regionData[x][y][z][c] = tmpData[tmpIndex + idx][c];
                        ++idx;
                    }
            idx = 0;
        }

        tmpIndex += eachNbin;
        ret.push_back(std::move(regionData));
    }
    return ret;
}




LoadDataResult loadOnflyData3D(const string& fnm) {
    LoadDataResult result;
    
    auto textData = readTextData(fnm);
    auto allData = readDataFile(fnm);
    
    // 解析头部 (基于MATLAB的textData索引，从0开始)
    if (textData.size() >= 2) {  // textData{2} -> index 1
        result.Lbox = extractNumbers(textData[1]);
    }
    if (textData.size() >= 3) {  // {3} -> index 2
        result.lowPos = extractNumbers(textData[2]);
    }
    if (textData.size() >= 4) {  // {4} -> index 3
        result.upPos = extractNumbers(textData[3]);
    }
    if (textData.size() >= 6) {  // {6} -> index 5
        auto nbin_nums = extractNumbers(textData[5]);
        result.nbin.resize(nbin_nums.size());
        for (size_t i = 0; i < nbin_nums.size(); ++i) {
            result.nbin[i] = static_cast<int>(nbin_nums[i]);
        }
    }
    
    if (textData.size() >= 8) {  // {8} -> index 7
        auto totNbin_nums = extractNumbers(textData[7]);
        if (!totNbin_nums.empty()) result.totNbin = static_cast<int>(totNbin_nums[0]);
    }
    if (textData.size() >= 9) {  // {9} -> index 8
        auto totFrames_nums = extractNumbers(textData[8]);
        if (!totFrames_nums.empty()) result.totFrames = static_cast<int>(totFrames_nums[0]);
    }
    
    result.nRegion = result.nbin.empty() ? 0 : result.nbin.size() / 3;
   
    
    // 生成xyzRange (linspace)
    for (size_t ii = 0; ii < result.nbin.size(); ++ii) {
        vector<double> range;
        double start;
        if (ii < result.lowPos.size()) {
            start = result.lowPos[ii];
        } else {
            std::cerr << "Error: lowPos size mismatch!" << std::endl;
            exit(1);
        }

        double end;
        if (ii < result.upPos.size()) {
            end = result.upPos[ii];
        } else {
            std::cerr << "Error: upPos size mismatch!" << std::endl;
            exit(1);
        }

        int points = result.nbin[ii];
        if (points > 1) {
            double step = (end - start) / (points - 1);
            for (int j = 0; j < points; ++j) {
                range.push_back(start + j * step);
            }
        } else {
            range.push_back(end);
        }
        result.xyzRange.push_back(range);
    }
    
    // 处理数据
    // if (result.nRegion != 1) {
    //     cerr << "Warning: Multiple regions not fully implemented. Assuming nRegion=1." << endl;
    // }
    int nx = (result.nRegion >= 1) ? result.nbin[0] : 1;
    int ny = (result.nRegion >= 1) ? result.nbin[1] : 1;
    int nz = (result.nRegion >= 1) ? result.nbin[2] : 1;
    int col = allData.empty() ? 0 : static_cast<int>(allData[0].size());
    
    result.ret.resize(result.totFrames);
    for (int ii = 0; ii < result.totFrames; ++ii) {
        int start_row = ii * result.totNbin;
        int end_row = min((ii + 1) * result.totNbin, static_cast<int>(allData.size()));
        
        vector<vector<double>> frameData;
        for (int j = start_row; j < end_row; ++j) {
            if (j < static_cast<int>(allData.size()) && static_cast<int>(allData[j].size()) == col) {
                frameData.push_back(allData[j]);
            }
        }
        
        result.ret[ii] = reshapeColumnMajor(frameData, result.nRegion, result.nbin, col);
    }
    result.nGroups = col / 4;  // Assuming 4 columns per group
    return result;
}

int main() {
    double V = 2.0;
    double scanrate = 0.0;
    
    string filename = std::string("./deal_data/onfly/onfly") + string(getenv("analyze_begin_case")) + "-" + 
                        string(getenv("analyze_end_case")) + ".dat";
    
    std::cout << "Loading data..." << endl;
    
    LoadDataResult data = loadOnflyData3D(filename);
    
    //添加对数组降维的支持
    if (data.ret.empty() || data.totFrames == 0) {
        cerr << "Error: No data loaded or unexpected format!" << endl;
        return -1;
    }
    
    // 使用ret (多帧4D数组)
    auto& retData = data.ret;
    remake_data charge_d(data.nRegion, data.nGroups), density_d(data.nRegion, data.nGroups);
    remake_data density(data.nRegion, 1), charge(data.nRegion, 1);

    int xindex = 0;
    int yindex = 0;//写成函数
    for (int region = 0; region < data.nRegion; ++region) {
        int nx = data.nbin[3 * region];       // 对应 nbinX
        int ny = data.nbin[3 * region + 1];   // 对应 nbinY
        int nz = data.nbin[3 * region + 2];   // 对应 nbinZ
        for (int n=0; n < data.nGroups; n++){
            for(int t=0; t<data.totFrames; t++){
                auto & tmp_charge_d = charge_d.get_region(region)[n].get_Matrix(nz, data.totFrames);
                auto & tmp_density_d = density_d.get_region(region)[n].get_Matrix(nz, data.totFrames);
                for(int z=0; z<nz;z++){
                    tmp_charge_d[z][t] = retData[t][region][xindex][yindex][z][n*4+2];
                    tmp_density_d[z][t] = retData[t][region][xindex][yindex][z][n*4+1];
                }
            }
            //
            if(n == 0) {
                density.get_region(region)[0] = density_d.get_region(region)[n];
                charge.get_region(region)[0] = charge_d.get_region(region)[n];
            }
            else{
                density.get_region(region)[0] += density_d.get_region(region)[n];
                charge.get_region(region)[0] += charge.get_region(region)[n];
            }
        }
    }

   //savedata();
    //几乎不可能储存成matlab的格式，至少以我的水平做不到。那这个的意义在哪里？

    return 0;
}