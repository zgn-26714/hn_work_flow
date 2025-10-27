#include "./include/matlab_to_C.h"
#include "./include/matrix.h"

using namespace std;
int main(int argc, char* argv[]) {
    string filename = argv[1];
    string begin_cas = getenv("analyze_begin_case");
    int len = begin_cas.size();
    string end_cas = getenv("analyze_end_case");
    itp::Matrix<double> data;
    vector<string> head;
    ifstream fid(filename);
    string line;
    while(getline(fid, line)){
        char tmp = line[0];
        if(tmp == '#' || tmp == '@')
            head.push_back(line);
        else
            break;
    }
    fid.close();

    for(int i=stoi(begin_cas); i<=stoi(end_cas); i++){
        itp::Matrix<double> tmp = matlab::xvgread(filename);
        if(i == stoi(begin_cas))
            data = tmp;
        else
         data += tmp;
    }
    for (int i=0; i<data.size(); i++){
        for (int j=0; j<data[0].size(); j++){
            data(i,j) /= (stoi(end_cas) - stoi(begin_cas) + 1);
        }
        cout << endl;
    }
    string out_dat = "ave_" + filename.substr(len,filename.size()-len-4) + begin_cas + "-" + end_cas + ".dat";
    matlab::xvgsave(out_dat, head, data.get());
}