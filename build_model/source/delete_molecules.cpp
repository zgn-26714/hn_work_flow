#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;

struct MolRecord {
    string mol_name;
    int old_count;
    int new_count;
    size_t line_idx;
    string indent;
};

struct PackmolEntry {
    string mol;
    int count;
};

struct AtomLine {
    string line;
};

struct MoleculeBlock {
    string resnr;
    string resname;
    vector<AtomLine> atoms;
};

// 去掉字符串首尾的空白字符，便于后续统一解析文本行。
static inline string trim(const string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// 按空白字符切分字符串，适合解析 packmol/top/gro 中的字段。
static inline vector<string> split(const string& s) {
    vector<string> tokens;
    istringstream iss(s);
    string tok;
    while (iss >> tok) tokens.push_back(tok);
    return tokens;
}

// 判断目标字符串是否存在于字符串数组中。
static inline bool contains(const vector<string>& v, const string& key) {
    for (const auto& x : v) {
        if (x == key) return true;
    }
    return false;
}

// 根据缩放因子计算每种目标分子的保留数量，并处理 FIXED_MOL 指定的同比例约束。
static void calc_factor(vector<MolRecord>& records, double factor) {
    string fixed_mol_env = getenv("FIXED_MOL") ? getenv("FIXED_MOL") : "";
    string force_scale_str = getenv("IS_FORCE_SCALE") ? getenv("IS_FORCE_SCALE") : "0";
    int is_force_scale_env = stoi(force_scale_str);

    unordered_set<string> fixed_set;
    {
        stringstream ss(fixed_mol_env);
        string token;
        while (ss >> token) {
            if (!token.empty()) fixed_set.insert(token);
        }
    }

    vector<MolRecord*> fixed_records;
    for (auto& rec : records) {
        if (fixed_set.count(rec.mol_name)) {
            fixed_records.push_back(&rec);
        }
    }

    if (!fixed_records.empty()) {
        sort(fixed_records.begin(), fixed_records.end(),
             [](const MolRecord* a, const MolRecord* b) {
                 return a->old_count < b->old_count;
             });

        int old_min = fixed_records.front()->old_count;
        int base = old_min;
        int step = 1;
        for (size_t i = 1; i < fixed_records.size(); ++i) {
            int oi = fixed_records[i]->old_count;
            int g = gcd(base, oi);
            int need = base / g;
            step = lcm(step, need);
        }

        int new_min;
        if (is_force_scale_env) {
            new_min = static_cast<int>(ceil(old_min * factor));
        } else {
            new_min = static_cast<int>(round(old_min * factor));
        }

        if (new_min % step != 0) {
            new_min = ((new_min / step) + 1) * step;
        }
        fixed_records.front()->new_count = new_min;

        for (size_t i = 1; i < fixed_records.size(); ++i) {
            fixed_records[i]->new_count =
                new_min * fixed_records[i]->old_count / old_min;
        }
    }

    for (auto& rec : records) {
        if (!fixed_set.count(rec.mol_name)) {
            if (is_force_scale_env) {
                rec.new_count = static_cast<int>(ceil(rec.old_count * factor));
            } else {
                rec.new_count = static_cast<int>(round(rec.old_count * factor));
            }
        }
    }
}

// 按 packmol 文件中的 structure/number 出现顺序读取分子及其数量。
static vector<PackmolEntry> read_packmol_sequence(const string& packmol_file) {
    ifstream fin(packmol_file);
    if (!fin) {
        throw runtime_error("Packmol file not found: " + packmol_file);
    }

    vector<PackmolEntry> seq;
    string line;
    string current_mol;
    bool in_structure = false;

    while (getline(fin, line)) {
        string s = trim(line);
        if (s.empty() || s[0] == '#') continue;

        auto parts = split(s);
        if (parts.size() >= 2 && parts[0] == "structure") {
            string pdb = parts[1];
            size_t pos = pdb.find_last_of('.');
            current_mol = (pos == string::npos) ? pdb : pdb.substr(0, pos);
            in_structure = true;
            continue;
        }

        if (in_structure && parts.size() == 2 && parts[0] == "number") {
            seq.push_back({current_mol, stoi(parts[1])});
            continue;
        }

        if (in_structure && parts.size() >= 2 && parts[0] == "end" && parts[1] == "structure") {
            in_structure = false;
            current_mol.clear();
        }
    }

    if (seq.empty()) {
        throw runtime_error("No structure/number blocks were found in " + packmol_file);
    }
    return seq;
}

// 修改 packmol 输入文件中目标分子的 number，并汇总修改后的总数量。
static void modify_packmol(const string& inp_file,
                           const vector<string>& target_mols,
                           double factor,
                           map<string, int>& total_target_counts) {
    ifstream fin(inp_file);
    if (!fin) {
        throw runtime_error("Input packmol file not found: " + inp_file);
    }

    vector<string> lines;
    string line;
    while (getline(fin, line)) {
        lines.push_back(line + "\n");
    }
    fin.close();

    vector<MolRecord> records;
    bool in_structure = false;
    string current_pdb;
    regex structure_re(R"(^\s*structure\s+([^#\s]+\.pdb))", regex::icase);

    for (size_t i = 0; i < lines.size(); ++i) {
        const string& raw = lines[i];
        string stripped = trim(raw);
        smatch match;

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

        if (in_structure && stripped.find("number") == 0) {
            auto parts = split(stripped);
            if (parts.size() < 2) continue;

            int old_count = 0;
            try {
                old_count = stoi(parts[1]);
            } catch (...) {
                continue;
            }

            string mol = regex_replace(current_pdb, regex(R"(\.pdb$)", regex::icase), "");
            if (!contains(target_mols, mol)) continue;

            size_t pos = raw.find("number");
            string indent = (pos != string::npos) ? raw.substr(0, pos) : "";
            records.push_back({mol, old_count, old_count, i, indent});
        }
    }

    if (records.empty()) {
        throw runtime_error(
            "No target molecules were matched for deletion. "
            "Please check MOL_name and packmol_.inp."
        );
    }

    calc_factor(records, factor);

    total_target_counts.clear();
    for (const auto& r : records) {
        total_target_counts[r.mol_name] += r.new_count;
        cout << "[packmol] " << r.mol_name << ": "
             << r.old_count << " -> " << r.new_count << " (delete branch)\n";
    }

    for (const auto& r : records) {
        ostringstream oss;
        oss << r.indent << "number " << r.new_count << "\n";
        lines[r.line_idx] = oss.str();
    }

    ofstream fout(inp_file);
    if (!fout) {
        throw runtime_error("Cannot write packmol file: " + inp_file);
    }
    for (const auto& l : lines) {
        fout << l;
    }
}

// 按新的 packmol 顺序和数量重写 topology 中 [ molecules ] 段。
static void modify_topology_from_packmol_sequence(
    const string& top_file,
    vector<PackmolEntry> packmol_seq) {
    ifstream fin(top_file);
    if (!fin) {
        throw runtime_error("Topology file not found: " + top_file);
    }

    vector<string> lines;
    string line;
    while (getline(fin, line)) {
        lines.push_back(line + "\n");
    }
    fin.close();

    vector<string> output;
    bool in_molecules = false;

    for (const auto& raw : lines) {
        string stripped = trim(raw);

        if (stripped == "[ molecules ]") {
            in_molecules = true;
            output.push_back(raw);
            continue;
        }

        if (in_molecules && !stripped.empty() && stripped[0] == '[' &&
            stripped.find("molecules") == string::npos) {
            in_molecules = false;
        }

        if (!in_molecules) {
            output.push_back(raw);
            continue;
        }

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
        for (auto it = packmol_seq.begin(); it != packmol_seq.end(); ++it) {
            if (it->mol == mol_top) {
                ostringstream oss;
                oss << "  " << left << setw(8) << mol_top
                    << right << setw(8) << it->count << "\n";
                output.push_back(oss.str());
                packmol_seq.erase(it);
                replaced = true;
                break;
            }
        }

        if (!replaced) {
            output.push_back(raw);
        }
    }

    ofstream fout(top_file);
    if (!fout) {
        throw runtime_error("Cannot write topology file: " + top_file);
    }
    for (const auto& l : output) {
        fout << l;
    }
}

// 读取 GRO 文件，并按分子块把原子记录重新组织起来。
static vector<MoleculeBlock> read_gro_molecules(const string& gro_file,
                                                string& header,
                                                string& box_line) {
    ifstream fin(gro_file);
    if (!fin) {
        throw runtime_error("GRO file not found: " + gro_file);
    }

    vector<MoleculeBlock> molecules;
    string line;
    if (!getline(fin, header)) {
        throw runtime_error("Invalid GRO file header: " + gro_file);
    }
    if (!getline(fin, line)) {
        throw runtime_error("Invalid GRO atom-count line: " + gro_file);
    }

    int atom_count = stoi(trim(line));
    MoleculeBlock current;
    bool has_current = false;

    for (int i = 0; i < atom_count; ++i) {
        if (!getline(fin, line)) {
            throw runtime_error("Unexpected EOF while reading GRO atoms: " + gro_file);
        }
        if (line.size() < 10) {
            throw runtime_error("Malformed GRO atom line: " + line);
        }

        string resnr = trim(line.substr(0, 5));
        string resname = trim(line.substr(5, 5));
        if (!has_current || current.resnr != resnr || current.resname != resname) {
            if (has_current) molecules.push_back(current);
            current = MoleculeBlock{};
            current.resnr = resnr;
            current.resname = resname;
            has_current = true;
        }
        current.atoms.push_back({line});
    }

    if (has_current) molecules.push_back(current);
    if (!getline(fin, box_line)) {
        throw runtime_error("Missing GRO box line: " + gro_file);
    }
    return molecules;
}

// 根据当前数量和目标数量，均匀挑选需要保留的分子下标。
static set<int> build_keep_indices(int current_count, int target_count) {
    set<int> keep;
    if (target_count < 0 || target_count > current_count) {
        throw runtime_error("Invalid target molecule count.");
    }
    if (target_count == current_count) {
        for (int i = 0; i < current_count; ++i) keep.insert(i);
        return keep;
    }
    if (target_count == 0) {
        return keep;
    }

    for (int k = 0; k < target_count; ++k) {
        int idx = static_cast<int>((static_cast<long long>(k) * current_count) / target_count);
        if (idx >= current_count) idx = current_count - 1;
        while (keep.count(idx) && idx < current_count - 1) {
            ++idx;
        }
        while (keep.count(idx) && idx > 0) {
            --idx;
        }
        keep.insert(idx);
    }
    return keep;
}

// 将筛选后的分子重新写回 GRO 文件，并更新原子总数。
static void write_gro(const string& gro_file,
                      const string& header,
                      const string& box_line,
                      const vector<MoleculeBlock>& molecules) {
    ofstream fout(gro_file);
    if (!fout) {
        throw runtime_error("Cannot write GRO file: " + gro_file);
    }

    int atom_count = 0;
    for (const auto& mol : molecules) {
        atom_count += static_cast<int>(mol.atoms.size());
    }

    fout << header << "\n";
    fout << atom_count << "\n";
    for (const auto& mol : molecules) {
        for (const auto& atom : mol.atoms) {
            fout << atom.line << "\n";
        }
    }
    fout << box_line << "\n";
}

// 按目标数量从 GRO 中删除指定分子，同时保留未参与删除的其他分子。
static void delete_from_gro(const string& gro_file,
                            const map<string, int>& target_totals,
                            const vector<string>& target_mols) {
    string header, box_line;
    auto molecules = read_gro_molecules(gro_file, header, box_line);

    map<string, vector<int>> indices_by_mol;
    for (int i = 0; i < static_cast<int>(molecules.size()); ++i) {
        if (contains(target_mols, molecules[i].resname)) {
            indices_by_mol[molecules[i].resname].push_back(i);
        }
    }

    set<int> keep_molecule_indices;
    for (int i = 0; i < static_cast<int>(molecules.size()); ++i) {
        if (!contains(target_mols, molecules[i].resname)) {
            keep_molecule_indices.insert(i);
        }
    }

    for (const auto& mol_name : target_mols) {
        auto idx_it = indices_by_mol.find(mol_name);
        if (idx_it == indices_by_mol.end()) continue;

        int current = static_cast<int>(idx_it->second.size());
        int desired = current;
        auto target_it = target_totals.find(mol_name);
        if (target_it != target_totals.end()) {
            desired = target_it->second;
        }
        if (desired > current) {
            throw runtime_error("Target count of " + mol_name +
                                " is larger than the current count in GRO.");
        }

        auto keep_local = build_keep_indices(current, desired);
        for (int local_idx : keep_local) {
            keep_molecule_indices.insert(idx_it->second[local_idx]);
        }

        cout << "[delete] " << mol_name << ": " << current << " -> " << desired << "\n";
    }

    vector<MoleculeBlock> kept;
    kept.reserve(keep_molecule_indices.size());
    for (int i = 0; i < static_cast<int>(molecules.size()); ++i) {
        if (keep_molecule_indices.count(i)) {
            kept.push_back(molecules[i]);
        }
    }

    write_gro(gro_file, header, box_line, kept);
}

int main(int argc, char* argv[]) {
    // 检查命令行参数个数是否正确。
    if (argc != 5) {
        // 参数数量不对时输出程序用法说明。
        cerr << "Usage: " << argv[0] << " SCALE_FACTOR INPUT_GRO TOPOLOGY_TOP PACKMOL_INP\n";
        // 参数错误时返回非零状态码。
        return 1;
    }

    // 读取第一个参数作为缩放因子，用于决定目标分子保留比例。
    double factor = atof(argv[1]);
    // 缩放因子不能为负数，否则直接报错退出。
    if (factor < 0.0) {
        // 提示用户 SCALE_FACTOR 的合法范围。
        cerr << "Error: SCALE_FACTOR must be non-negative\n";
        // 非法输入时返回错误码。
        return 1;
    }

    // 读取输入的 GRO 文件路径。
    const string gro_file = argv[2];
    // 读取输入的 topology 文件路径。
    const string top_file = argv[3];
    // 读取输入的 packmol 文件路径。
    const string packmol_file = argv[4];

    // 从环境变量中获取需要删除或缩放的目标分子名称列表。
    const char* mol_env = getenv("MOL_name");
    // 如果没有设置 MOL_name，则无法知道要处理哪些分子。
    if (!mol_env) {
        // 输出缺少环境变量的错误信息。
        cerr << "Error: environment variable MOL_name must be set\n";
        // 环境变量缺失时返回错误码。
        return 1;
    }

    // 将环境变量中的分子名按空白切分成数组。
    vector<string> target_mols = split(mol_env);
    // 若解析后为空，说明 MOL_name 没有提供有效分子名。
    if (target_mols.empty()) {
        // 输出空分子列表的错误信息。
        cerr << "Error: MOL_name is empty\n";
        // 输入无效时返回错误码。
        return 1;
    }

    // 捕获整个处理流程中的异常，统一输出错误信息。
    try {
        // 保存 packmol 修改后每种目标分子的总数量，供后续 GRO 删除步骤使用。
        map<string, int> total_target_counts;
        // 先修改 packmol 文件中的目标分子数量。
        modify_packmol(packmol_file, target_mols, factor, total_target_counts);
        // 重新读取 packmol 中 structure/number 的顺序和数量。
        auto packmol_seq = read_packmol_sequence(packmol_file);
        // 按最新的 packmol 结果同步更新 topology 的 [ molecules ] 段。
        modify_topology_from_packmol_sequence(top_file, packmol_seq);
        // 按目标数量从 GRO 文件中删除对应分子。
        delete_from_gro(gro_file, total_target_counts, target_mols);
        // 输出 GRO 文件已更新的提示信息。
        cout << "Updated GRO written to: " << gro_file << "\n";
        // 输出 topology 文件已更新的提示信息。
        cout << "Updated topology written to: " << top_file << "\n";
        // 输出 packmol 文件已更新的提示信息。
        cout << "Updated packmol written to: " << packmol_file << "\n";
    // 捕获标准异常并输出具体错误内容。
    } catch (const exception& e) {
        // 输出异常描述，便于定位失败原因。
        cerr << "Error: " << e.what() << "\n";
        // 异常情况下返回错误码。
        return 1;
    }

    // 所有步骤成功完成后返回 0。
    return 0;
}
