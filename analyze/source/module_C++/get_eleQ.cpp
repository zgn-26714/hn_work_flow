#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;
using real = double;

int main(int argc, char *argv[]){
    string begin = getenv("analyze_begin_case");
    string end = getenv("analyze_end_case");
    string line;
    string molecule = getenv("analyze_mol");
    vector<real> charge;
    float qout, md_dt = stof(getenv("DT")); // Default value for md_dt
    string begin_dir = "./case" + begin + "/CPM_ControlFile.dat";
    ifstream CPM_Control(begin_dir);
    if (!CPM_Control.is_open()) {
        cerr << "not exist " << begin_dir << endl;
        return 1;
    }
    while (getline(CPM_Control, line))
    {
        if (line[0] == 'f'){
            qout = stoi(line.substr(line.find('=') + 1,line.find('#'))) * md_dt;
        }
    }
    CPM_Control.close();

    for (int i = stoi(begin); i <= stoi(end); i++){
        string CPMfile = "./case" + to_string(i) + "/CPM_electrodeCharge.dat";
        cout<<"\rcase is "<<i<<flush;
        ifstream file(CPMfile);
        if (!file.is_open()) {
            cerr << "Error opening file: " << CPMfile << endl;
            continue;
        }
        else{
            int count = 0;
            while(getline(file, line)){
                if(line[0] == '#') continue; // Skip comment lines
                while (line[0] == ' ') line.erase(0, 1); // Remove leading spaces
                real tmp = stod(line.substr(0, line.find(' '))); // Extract the first double value
                // if ( tmp < 0) tmp = -tmp; 
                if(i == stoi(begin)){
                    charge.emplace_back(tmp);
                } 
                else {
                    charge[count] += tmp;
                }
                count++;
            }
        }
    }
    cout<<endl;
    for (int i = 0; i < charge.size(); i++){
        charge[i] /= (stoi(end) - stoi(begin) + 1);
    }
    
    string outfile = "./deal_data/eleQ/" + molecule + "_eleCharge" + string(getenv("analyze_V")) + "V" +
                                string(getenv("analyze_tau")) + "ps_" + begin + "-" + end + ".dat";
    ofstream output(outfile);
    
    if (!output) {
        cerr << "Error: Cannot open file " << outfile << endl;
        return 1;
    }
    
    output << "# ElecCharge averaged from case " << begin << " to " << end << endl;
    output << "# V = " << string(getenv("analyze_V")) << " V" << endl;
    output << "# tau = " << string(getenv("analyze_tau")) << " ps" << endl;
    output << "# MD dt = " << string(getenv("DT")) << " ps" << endl;
    output << "# Unit: e" << endl;
    output << "# qout = " << qout << endl;
    for (int i = 0; i < charge.size(); i++){
        output << -charge[i] << endl;
    }
    output.close();
    cout << "Output written to " << outfile << endl;
    return 0;
}
