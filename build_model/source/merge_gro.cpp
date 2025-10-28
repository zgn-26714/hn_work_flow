#include <string>
#include <fstream>
#include <vector>
#include <sstream>
using namespace std;

const string ele_name = getenv("ele_name");
int main(int argc, char* argv[]) {
    string gro1 = argv[1];
    string gro2 = argv[2];
    ifstream file1(gro1);
    ifstream file2(gro2);
    vector<string> pre_head;
    vector<string> merge_ele;
    vector<string> merge_lines;
    string pre_end;
    string line;
    int count = 0;
    while (getline(file1, line)) {
        if(count < 2){
            pre_head.push_back(line);
            count++;
            continue;
        }
        std::istringstream iss(line);
        string first;
        iss>>first;
        string mol_name;
        for(char i : first){
            if (i<'0' || i>'9') 
                mol_name += i; 
        }
        if(mol_name.empty()){
            pre_end = line;
            continue;
        }
        if(ele_name.find(mol_name) == string::npos)
            merge_lines.push_back(line);
        else
            merge_ele.push_back(line);
    }
    file1.close();

    while (getline(file2, line)) {
        if(count < 2){
            pre_head.push_back(line);
            count++;
            continue;
        }
        std::istringstream iss(line);
        string first;
        iss>>first;
        string mol_name;
        for(char i : first){
            if (i<'0' || i>'9') 
                mol_name += i; 
        }
        if(mol_name.empty()){
            pre_end = line;
            continue;
        }
        if(ele_name.find(mol_name) == string::npos)
            merge_lines.push_back(line);
        else
            merge_ele.push_back(line);
    }
    file2.close();


    int num = stoi(pre_head[1] + pre_head[3]);

    string outfile_name = gro1.substr(0, gro1.size()-4) + "_merge.gro";
    ofstream outfile(outfile_name);
    outfile << pre_head[0] << endl;
    outfile << num << endl;
    for(string l : merge_ele)
        outfile << l << endl;
    for(string l : merge_lines)
        outfile << l << endl;
    outfile << pre_end[0];



}