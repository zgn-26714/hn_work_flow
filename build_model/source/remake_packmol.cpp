#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <unordered_set>
#include <numeric>

using namespace std;

struct MolRecord {
    string mol_name;     // 分子名（去掉 .pdb）
    int old_count;       // 原始 number
    int new_count;       // 计算后的 number
    size_t line_idx;     // number 行在 lines 中的位置
    string indent;       // 原始缩进
};

struct PackmolEntry {
    string mol;   // 由 structure xxx.pdb 得到
    int count;    // number N
};


/* ================= Utility ================= */
void calc_factor(vector<MolRecord>& records, double factor) {
    const char* fixed_mol_env = getenv("FIXED_MOL");
    int is_force_scale_env = stoi(string(getenv("IS_FORCE_SCALE")));
    if (!fixed_mol_env || !is_force_scale_env)
        throw std::runtime_error("Environment variables not set");
    
    std::string fixed_mol(fixed_mol_env);
    std::unordered_set<std::string> fixed_set;
    {
        std::stringstream ss(fixed_mol);
        std::string token;

        while (ss >> token) {
            if (!token.empty())
                fixed_set.insert(token);
        }
    }
     /* ---------- 2. 收集 FIXED_MOL 对应的 records ---------- */
    std::vector<MolRecord*> fixed_records;
    for (auto& rec : records) {
        if (fixed_set.count(rec.mol_name)) {
            fixed_records.push_back(&rec);
        }
    }

    /* ---------- 3. FIXED_MOL 内部分子：保持比例 ---------- */
    if (!fixed_records.empty()) {

        // 按 old_count 从小到大排序
        std::sort(fixed_records.begin(), fixed_records.end(),
                  [](const MolRecord* a, const MolRecord* b) {
                      return a->old_count < b->old_count;
                  });

        int old_min = fixed_records.front()->old_count;
        // 其他 FIXED_MOL 分子：按倍数关系
        int base = old_min;

        int step = 1;
        for (size_t i = 1; i < fixed_records.size(); ++i) {
            int oi = fixed_records[i]->old_count;
            int g  = std::gcd(base, oi);
            int need = base / g;
            step = std::lcm(step, need);
        }

        int new_min;
        if (is_force_scale_env){
            new_min = static_cast<int>(std::ceil(old_min * factor));
        }
        else{
            new_min = static_cast<int>(std::round(old_min * factor));
        }

        // 向上调整到满足比例的最小值
        if (new_min % step != 0) {
            new_min = ((new_min / step) + 1) * step;
        }

        fixed_records.front()->new_count = new_min;

        for (size_t i = 1; i < fixed_records.size(); ++i) {
            fixed_records[i]->new_count =
                new_min * fixed_records[i]->old_count / old_min;
        }
    }

    /* ---------- 4. 非 FIXED_MOL 分子：简单乘 factor ---------- */
    for (auto& rec : records) {
        if (!fixed_set.count(rec.mol_name)) {
            if (is_force_scale_env)
                rec.new_count =
                    static_cast<int>(std::ceil(rec.old_count * factor));
            else    
                rec.new_count =
                    static_cast<int>(std::round(rec.old_count * factor));
        }
    }

}




static inline string trim(const string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static inline vector<string> split(const string& s) {
    vector<string> tokens;
    istringstream iss(s);
    string tok;
    while (iss >> tok) tokens.push_back(tok);
    return tokens;
}

static inline bool contains(const vector<string>& v, const string& key) {
    for (const auto& x : v)
        if (x == key) return true;
    return false;
}

vector<PackmolEntry>
read_packmol_sequence(const string& packmol_file)
{
    ifstream fin(packmol_file);
    if (!fin)
        throw runtime_error("Packmol file not found: " + packmol_file);

    vector<PackmolEntry> seq;

    string line;
    string current_mol;
    bool in_structure = false;

    while (getline(fin, line)) {
        string s = trim(line);
        if (s.empty() || s[0] == '#')
            continue;

        auto parts = split(s);

        /* structure xxx.pdb */
        if (parts.size() >= 2 && parts[0] == "structure") {
            string pdb = parts[1];
            size_t pos = pdb.find_last_of('.');
            current_mol = (pos == string::npos)
                          ? pdb
                          : pdb.substr(0, pos);
            in_structure = true;
            continue;
        }

        /* number N */
        if (in_structure && parts.size() == 2 && parts[0] == "number") {
            int n = stoi(parts[1]);
            seq.push_back({current_mol, n});
            continue;
        }

        /* end structure */
        if (in_structure &&
            parts.size() >= 2 &&
            parts[0] == "end" &&
            parts[1] == "structure") {
            in_structure = false;
            current_mol.clear();
        }
    }

    fin.close();

    if (seq.empty()){
        cerr << "[packmol] Warning: no structures found\n";
        exit(0);
    }

    return seq;
}


/* ================= Packmol ================= */

void modify_packmol(const string& inp_file,
                    const vector<string>& target_mols,
                    double factor)
{
    ifstream fin(inp_file);
    if (!fin)
        throw runtime_error("Input packmol file not found: " + inp_file);

    /* ---------- 读文件 ---------- */
    vector<string> lines;
    string line;
    while (getline(fin, line))
        lines.push_back(line + "\n");
    fin.close();

    /* ---------- 第一阶段：解析 ---------- */
    vector<MolRecord> records;

    bool in_structure = false;
    string current_pdb;

    regex structure_re(R"(^\s*structure\s+([^#\s]+\.pdb))",
                       regex::icase);

    for (size_t i = 0; i < lines.size(); ++i) {
        const string& raw = lines[i];
        string stripped = trim(raw);
        smatch match;

        /* structure xxx.pdb */
        if (regex_search(stripped, match, structure_re)) {
            in_structure = true;
            current_pdb = match[1];
            continue;
        }

        if (in_structure && stripped.find("end structure") == 0) {
            in_structure = false;
            current_pdb.clear();
            continue;
        }

        /* number 行 */
        if (in_structure && stripped.find("number") == 0) {
            auto parts = split(stripped);
            if (parts.size() < 2)
                continue;

            int old_count;
            try {
                old_count = stoi(parts[1]);
            } catch (...) {
                continue;
            }

            string mol = regex_replace(current_pdb,
                                       regex(R"(\.pdb$)", regex::icase),
                                       "");

            if (!contains(target_mols, mol))
                continue;

            size_t pos = raw.find("number");
            string indent = (pos != string::npos) ? raw.substr(0, pos) : "";

            records.push_back({
                mol,
                old_count,
                old_count,   // new_count 先占位
                i,
                indent
            });
        }
    }

    /* ---------- 第二阶段：计算 ---------- */
    calc_factor(records, factor);
    for (const auto& r : records) {
        std::cout << "[packmol] " << r.mol_name << ": "
             << r.old_count << " -> "
             << r.new_count << " (x"
             << fixed << setprecision(3)
             << (static_cast<double>(r.new_count) / r.old_count)
             << ")\n";
    }
    /* ---------- 第三阶段：回写 ---------- */
    for (const auto& r : records) {
        ostringstream oss;
        oss << r.indent << "number " << r.new_count << "\n";
        lines[r.line_idx] = oss.str();
    }

    ofstream fout(inp_file);
    for (const auto& l : lines)
        fout << l;
    fout.close();

    std::cout << "Modified packmol input written to: "
         << inp_file << endl;
}


/* ================= Topology ================= */

void modify_topology_from_packmol_sequence(
    const string& top_file,
    vector<PackmolEntry> packmol_seq)
{
    ifstream fin(top_file);
    if (!fin)
        throw runtime_error("Topology file not found: " + top_file);

    vector<string> lines;
    string line;
    while (getline(fin, line))
        lines.push_back(line + "\n");
    fin.close();

    vector<string> output;
    bool in_molecules = false;

    for (const auto& raw : lines) {
        string stripped = trim(raw);

        /* enter [ molecules ] */
        if (stripped == "[ molecules ]") {
            in_molecules = true;
            output.push_back(raw);
            continue;
        }

        /* leave [ molecules ] */
        if (in_molecules &&
            !stripped.empty() &&
            stripped[0] == '[' &&
            stripped.find("molecules") == string::npos) {
            in_molecules = false;
        }

        if (!in_molecules) {
            output.push_back(raw);
            continue;
        }

        /* comments / empty */
        if (stripped.empty() || stripped[0] == ';') {
            output.push_back(raw);
            continue;
        }

        auto parts = split(stripped);
        if (parts.size() < 2) {
            output.push_back(raw);
            continue;
        }

        string mol_top = parts[0];
        bool replaced = false;

        /* 在 packmol 顺序表中查找同名分子 */
        for (auto it = packmol_seq.begin(); it != packmol_seq.end(); ++it) {
            if (it->mol == mol_top) {
                int new_count = it->count;
                if (new_count < 0)
                    throw runtime_error("Negative molecule count for " + mol_top);

                ostringstream oss;
                oss << "  " << left << setw(8) << mol_top
                    << right << setw(8) << new_count << "\n";
                output.push_back(oss.str());

                cout << "[top] " << mol_top
                     << " <- " << new_count
                     << " (from packmol)\n";

                /* 消费掉这个 packmol 记录，避免重复使用 */
                packmol_seq.erase(it);
                replaced = true;
                break;
            }
        }

        if (!replaced) {
            /* packmol 中没找到同名分子，原样保留 */
            output.push_back(raw);
        }
    }

    ofstream fout(top_file);
    for (const auto& l : output)
        fout << l;
    fout.close();

    cout << "Modified topology written to: "
         << top_file << endl;
}



/* ================= main ================= */

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " SCALE_FACTOR packmol_name\n";
        return 1;
    }

    double factor = atof(argv[1]);
    if (factor < 0.0) {
        cerr << "Error: SCALE_FACTOR must be non-negative\n";
        return 1;
    }

    const char* iter = getenv("iter");
    const char* top_env = getenv("TOP");
    const char* mol_env = getenv("MOL_name");

    if (!iter || !top_env || !mol_env) {
        cerr << "Error: environment variables iter / TOP / MOL_name must be set\n";
        return 1;
    }

    vector<string> target_mols = split(mol_env);

    std::cout << "\n"
         << " ---------------------------------------- iter=\"" << iter << "\" -----------------------------------------\n"
         << " [INFO] change molecule num...\n"
         << " Script    : " << argv[0] << "\n";

    try {
        modify_packmol(string(argv[2]) + ".inp",
                        target_mols, factor);

        auto packmol_seq = read_packmol_sequence(string(argv[2]) + ".inp");

        modify_topology_from_packmol_sequence(string(top_env) + ".top", packmol_seq);
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    std::cout << "✅ All files modified successfully!\n";
    return 0;
}
