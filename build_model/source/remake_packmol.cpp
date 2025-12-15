#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <cstdlib>
#include <cmath>
#include <iomanip>

using namespace std;

/* ================= Utility ================= */

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

/* ================= Packmol ================= */

void modify_packmol(const string& inp_file,
                    const vector<string>& target_mols,
                    double factor)
{
    struct MolCount {
        int origin;
        int new_count;
    };
    vector<MolCount> data_mol_nums;
    pair<int, float> origin_ratio;
    ifstream fin(inp_file);
    if (!fin)
        throw runtime_error("Input packmol file not found: " + inp_file);

    vector<string> lines;
    string line;
    while (getline(fin, line))
        lines.push_back(line + "\n");
    fin.close();

    vector<string> output;
    bool in_structure = false;
    string current_pdb;

    regex structure_re(R"(^\s*structure\s+([^#\s]+\.pdb))",
                       regex::icase);

    for (const auto& raw : lines) {
        string stripped = trim(raw);
        smatch match;

        /* structure xxx.pdb */
        if (regex_search(stripped, match, structure_re)) {
            in_structure = true;
            current_pdb = match[1];
            output.push_back(raw);
            continue;
        }

        if (in_structure && stripped.find("end structure") == 0) {
            in_structure = false;
            current_pdb.clear();
            output.push_back(raw);
            continue;
        }

        if (in_structure && stripped.find("number") == 0) {
            auto parts = split(stripped);
            if (parts.size() < 2) {
                cerr << "[packmol] Warning: invalid number line: "
                     << stripped << endl;
                output.push_back(raw);
                continue;
            }

            int old_count;
            try {
                old_count = stoi(parts[1]);
            } catch (...) {
                cerr << "[packmol] Warning: invalid number value: "
                     << stripped << endl;
                output.push_back(raw);
                continue;
            }

            string mol = regex_replace(current_pdb,
                                       regex(R"(\.pdb$)", regex::icase),
                                       "");

            if (contains(target_mols, mol)) {
                int new_count = int(round(old_count * factor));

                if (new_count < 0)
                    throw runtime_error("Negative molecule count for " + mol);

                size_t pos = raw.find("number");
                string indent = (pos != string::npos) ? raw.substr(0, pos) : "";

                ostringstream oss;
                oss << indent << "number " << new_count << "\n";
                output.push_back(oss.str());

                cout << "[packmol] " << mol << ": "
                     << old_count << " -> "
                     << new_count << " (x"
                     << fixed << setprecision(3) << factor << ")\n";
            } else {
                output.push_back(raw);
            }
        }
        else {
            output.push_back(raw);
        }
    }

    ofstream fout(inp_file);
    for (const auto& l : output)
        fout << l;
    fout.close();

    cout << "Modified packmol input written to: "
         << inp_file << endl;
}

/* ================= Topology ================= */

void modify_topology(const string& top_file,
                     const vector<string>& target_mols,
                     double factor)
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

        if (stripped == "[ molecules ]") {
            in_molecules = true;
            output.push_back(raw);
            continue;
        }

        if (in_molecules &&
            stripped.size() > 0 &&
            stripped[0] == '[' &&
            stripped.find("molecules") == string::npos) {
            in_molecules = false;
        }

        if (in_molecules) {
            if (stripped.empty() || stripped[0] == ';') {
                output.push_back(raw);
                continue;
            }

            auto parts = split(stripped);
            if (parts.size() < 2) {
                output.push_back(raw);
                continue;
            }

            string mol = parts[0];
            int old_count;

            try {
                old_count = stoi(parts[1]);
            } catch (...) {
                cerr << "[top] Warning: invalid count: "
                     << stripped << endl;
                output.push_back(raw);
                continue;
            }

            if (contains(target_mols, mol)) {
                int new_count = int(round(old_count * factor));
                if (new_count < 0)
                    throw runtime_error("Negative molecule count for " + mol);

                ostringstream oss;
                oss << "  " << left << setw(8) << mol
                    << right << setw(8) << new_count << "\n";
                output.push_back(oss.str());

                cout << "[top] " << mol << ": "
                     << old_count << " -> "
                     << new_count << " (x"
                     << fixed << setprecision(3) << factor << ")\n";
            } else {
                output.push_back(raw);
            }
        }
        else {
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

    cout << "\n"
         << "===================================================================================================\n"
         << " ---------------------------------------- iter=\"" << iter << "\" -----------------------------------------\n"
         << " [INFO] Starting pre-equilibration and density analysis workflow\n"
         << " Script    : " << argv[0] << "\n"
         << "===================================================================================================\n";

    try {
        modify_packmol(string(argv[2]) + ".inp",
                        target_mols, factor);

        modify_topology(string(top_env) + ".top",
                        target_mols, factor);
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    cout << "✅ All files modified successfully!\n";
    return 0;
}
