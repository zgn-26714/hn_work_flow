#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/stat.h>
using namespace std;

// Global parameters
int _reserve_num = 0;                   // reserved vector capacity
int num = 100;                          // number of output frames
int T = 298;                            // target temperature
string input;                           // input file prefix
string outdir;                          // output directory
string outenergy;                       // output energy file
string command1 = "echo temperature | gmx energy -f ?.edr -o ? > test.log 2>&1";   // command to extract temperature from .edr
string command2 = "echo 0 | gmx trjconv -f ?.trr -s ?.tpr -b ? -e ? -o ? > test.log 2>&1  "; // command to extract frames
vector<double> temps;                   // temperature deviations
vector<double> times;                   // corresponding times

// Check if file exists
bool file_exists(const string& name) {
    ifstream f(name.c_str());
    return f.good();
}

// Print help message
void print_help() {
    cout << "Usage: hn -i input_prefix -t temperature -n num -o output_dir\n"
         << "Options:\n"
         << "  -i   Input file prefix (must include .tpr, .trr, .edr)\n"
         << "  -t   Target temperature (integer)\n"
         << "  -n   Number of frames to extract (default: 100)\n"
         << "  -o   Output directory\n"
         << "  -h   Print this help message\n"
         << endl;
}

int main(int argc, char* argv[]) {
    // If no arguments, show help
    if (argc == 1) {
        print_help();
        return 0;
    }

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) input = argv[++i];
        else if (arg == "-t" && i + 1 < argc) T = stoi(argv[++i]);
        else if (arg == "-n" && i + 1 < argc) num = stoi(argv[++i]);
        else if (arg == "-o" && i + 1 < argc) outdir = argv[++i];
        else if (arg == "-h") {
            print_help();
            return 0;
        }
    }

    // Basic checks
    if (input.empty() || outdir.empty()) {
        cerr << "Error: Missing required arguments!\n" << endl;
        print_help();
        return 1;
    }

    // Check required input files
    if (!file_exists(input + ".tpr") || !file_exists(input + ".trr") || !file_exists(input + ".edr")) {
        cerr << "Error: Missing input files (.tpr, .trr, .edr) with prefix: " << input << endl;
        return 1;
    }

    // Create output directory
#ifdef _WIN32
    system(("mkdir " + outdir).c_str());
#else
    system(("mkdir -p " + outdir).c_str());
#endif

    outenergy = outdir + "/temp.xvg";

    // Run gmx energy command to extract temperature
    string c1 = command1;
    c1.replace(c1.find('?'), 1, input);
    c1.replace(c1.find('?'), 1, outenergy);
    system(c1.c_str());

    // Reserve space for vectors
    if (_reserve_num == 0) _reserve_num = num * 100;
    temps.reserve(_reserve_num);
    times.reserve(_reserve_num);

    // Read temp.xvg
    ifstream readfile(outenergy);
    if (!readfile.is_open()) {
        cerr << "Error: Cannot open " << outenergy << endl;
        return 1;
    }

    string line;
    while (getline(readfile, line)) {
        if (line.empty()) continue;
        if (line[0] == '@' || line[0] == '#') continue;
        while (!line.empty() && line[0] == ' ') line = line.substr(1);
        size_t pos = line.find(' ');
        if (pos == string::npos) continue;
        times.push_back(stod(line.substr(0, pos)));
        line = line.substr(pos + 1);
        while (!line.empty() && line[0] == ' ') line = line.substr(1);
        temps.push_back(stod(line));
    }
    readfile.close();

    // Convert temps to deviation from T
    for (auto & temp : temps) temp = abs(temp - T);

    // Sort temps and times based on deviation
    for (int i = 1; i < (int)temps.size(); i++) {
        if (temps[i] >= temps[i - 1]) continue;
        auto index = lower_bound(temps.begin(), temps.begin() + i - 1, temps[i]) - temps.begin();
        temps.insert(temps.begin() + index, temps[i]);
        temps.erase(temps.begin() + i + 1);
        times.insert(times.begin() + index, times[i]);
        times.erase(times.begin() + i + 1);
    }

    // Write time-temperature pairs
    ofstream tandt(outdir + "/tandt.txt");
    if (!tandt.is_open()) {
        cerr << "Error: Cannot write tandt.txt" << endl;
        return 1;
    }
    for (int i = 0; i < (int)times.size(); i++)
        tandt << times[i] << '\t' << temps[i] << endl;
    tandt.close();

    // Prepare gmx trjconv command
    string c2base = command2;
    c2base.replace(c2base.find('?'), 1, input);
    c2base.replace(c2base.find('?'), 1, input);

    // Extract frames
    for (int i = 0; i < num && i < (int)times.size(); i++) {
        string getframes(c2base);
        string output = outdir + "/frame" + to_string(i + 1) + ".gro";
        getframes.replace(getframes.find('?'), 1 , to_string(times[i] - 0.001));
        getframes.replace(getframes.find('?'), 1 , to_string(times[i] + 0.001));
        getframes.replace(getframes.find('?'), 1 , output);
        system(getframes.c_str());
        std::cout<<"[frame"<<i<<"] success!"<<std::endl;
    }

    return 0;
}
