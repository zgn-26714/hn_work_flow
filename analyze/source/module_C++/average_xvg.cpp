#include "./include/matlab_to_C.h"
#include "./include/matrix.h"

using namespace std;

void checkFileExist(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: File " << filename << " does not exist." << endl;
        exit(1);
    }
    file.close();
}

void fill_head(const string& filename, vector<string>& head) {
    ifstream fid(filename);

    string str;
    while (getline(fid, str)) {
        if (str.empty()) continue;
        char tmp = str[0];
        if (tmp == '#' || tmp == '@') {
            head.push_back(str);
        } else {
            break;
        }
    }
    fid.close();
}


int main() {
    string analyze_cpp = getenv("analyze_cpp");
    string file_fold = "./deal_data/" + analyze_cpp + "/";
    string begin_cas = getenv("analyze_begin_case");
    string mol_name = getenv("analyze_mol");
    string end_cas = getenv("analyze_end_case");
    string is_skip_first_line = getenv("analyze_is_skip_first_line");
    itp::Matrix<double> data;
    vector<string> head;
    
    for(int i=stoi(begin_cas); i<=stoi(end_cas); i++){
        string filename = file_fold + to_string(i) + analyze_cpp + ".xvg";
        checkFileExist(filename);
        if(head.empty()){
            fill_head(filename, head);
        }
        itp::Matrix<double> tmp = matlab::xvgread(filename);
        if(i == stoi(begin_cas))
            data = tmp;
        else{
            if (data.size() != tmp.size() || data[0].size() != tmp[0].size() || data[1].size() != tmp[1].size()){
                cerr<<"Error: File "<<filename<<" has different data size!"<<endl;
                cerr<<"size of pre_Data: "<<data.size()<<", size of curr_Data: "<<tmp.size()<<endl;
                cerr<<"size of pre_Data first line: "<<data[0].size()<<", size of curr_Data first line: "<<tmp[0].size()<<endl;
                cerr<<"size of pre_Data other line: "<<data[1].size()<<", size of curr_Data other line: "<<tmp[1].size()<<endl;
                exit(1);
            }
             data += tmp;
        }
    }
    for (int i=0; i<data.size(); i++){
        if(i == 0 && stoi(is_skip_first_line))
            continue;
        for (int j=0; j<data[i].size(); j++){
            data(i,j) /= (stoi(end_cas) - stoi(begin_cas) + 1);
        }
    }
    string out_dat = file_fold + mol_name +"_ave_" + analyze_cpp + "_" + begin_cas + "-" + end_cas + ".dat";
    cout<<"saving to "<<out_dat<<endl;
    matlab::xvgsave(out_dat, head, data.get());
}