#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
using namespace std;

void changeMol_name(string& mol_lines);

int main(int argc, char* argv[]) {
    string gro1 = argv[1];
    string gro2 = argv[2];
    ifstream file1(gro1);
    ifstream file2(gro2);
    vector<string> pre_head;
    vector<string> all_name;
    vector<vector<string>> merge_lines;
    string pre_end;
    string line;
    
    vector<string> mol_lines;

    if (!file1.is_open() || !file2.is_open()) {
        cerr << "Error: Cannot open input files" << endl;
        return 1;
    }

    string now_name;
    int count = 0;
    while (getline(file1, line)) {
        
        if(count < 2){
            pre_head.push_back(line);
            count++;
            continue;
        }
        std::istringstream iss(line);
        string mol_name;
        iss>>mol_name;
        size_t pos = mol_name.find_first_not_of("0123456789");
        if (pos != string::npos) {
            mol_name = mol_name.substr(pos);
        }

        if(mol_name[0] == '.'){
            pre_end = line;
            continue;
        }
        //new mol stucture
        if(mol_name != now_name){
            auto iter = find(all_name.begin(),all_name.end(),now_name);
            //new mol
            if(iter == all_name.end()){
                if(!now_name.empty())
                    all_name.push_back(now_name);
                if(!mol_lines.empty()){
                    merge_lines.push_back(mol_lines);
                    mol_lines.clear();
                }
            }
            else{//have mol with same name
                int index = iter - all_name.begin();
                if(!mol_lines.empty()){
                    merge_lines[index].insert(merge_lines[index].end(), mol_lines.begin(), mol_lines.end());
                    mol_lines.clear();
                }
            }
            now_name = mol_name;
        }

            mol_lines.push_back(line);

    }
    file1.close();
    auto iter = find(all_name.begin(),all_name.end(),now_name);
    if(iter == all_name.end()){
        if(!merge_lines.empty()){
            all_name.push_back(now_name);
            merge_lines.push_back(mol_lines);
            mol_lines.clear();
        }
    }
    else{
        int index = iter - all_name.begin();
        if(!mol_lines.empty()){
            merge_lines[index].insert(merge_lines[index].end(), mol_lines.begin(), mol_lines.end());
            mol_lines.clear();
        }
    }
    std::cout<<"Finished reading "<< gro1 << std::endl;

    count = 0;
    while (getline(file2, line)) {
        if(count < 2){
            pre_head.push_back(line);
            count++;
            continue;
        }
        std::istringstream iss(line);
        string mol_name;
        iss>>mol_name;
        size_t pos = mol_name.find_first_not_of("0123456789");
        if (pos != string::npos) {
            mol_name = mol_name.substr(pos);
        }
        if(mol_name[0] == '.'){
            pre_end = line;
            continue;
        }
        auto iter = find(all_name.begin(),all_name.end(),mol_name);
        if(iter == all_name.end()){
            cerr << "Error: Molecule " << mol_name << " in " << gro2 << " not found in " << gro1 << endl;
            return 1;
        }
        else if(iter == all_name.begin()){
            changeMol_name(line);
            merge_lines[0].push_back(line);
        }
        else{
            int index = iter - all_name.begin();
            merge_lines[index].push_back(line);
        }
        
    }
    file2.close();
    

    int num = stoi(pre_head[1]) + stoi(pre_head[3]);
    string outfile_name = gro1.substr(0, gro1.size()-4) + "_merge.gro";
    ofstream outfile(outfile_name);
    outfile << pre_head[0] << endl;
    outfile << num << endl;
    count = 0;
    for(auto& mol_lines : merge_lines){
        for(string l : mol_lines)
            outfile << l << endl;
    }
    outfile << pre_end;
}

void changeMol_name(string& line){
    std::istringstream iss(line);
    string tmp;
    iss >> tmp;
    int num = line.find(tmp);
    line = line.substr(num + tmp.size(), line.size());
    string space(num, ' ');
    size_t pos = tmp.find_first_not_of("0123456789");
    int space_add = tmp.size() - pos;
    if (pos != string::npos) {
        tmp = tmp.substr(0, pos);
    }
    string space2;
    if (space_add > 2)
        space2 = string(space_add - 2, ' ');
    line = space + tmp + "ER" + space2 + line;
}