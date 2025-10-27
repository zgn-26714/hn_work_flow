#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <regex>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input_file> [output_file]" << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    std::string outputFile = (argc >= 3) ? argv[2] : (inputFile + ".processed");
    
    std::ifstream infile(inputFile);
    std::ofstream outfile(outputFile);
    
    if (!infile.is_open() || !outfile.is_open()) {
        std::cerr << "Error: Cannot open file" << std::endl;
        return 1;
    }
    
    std::unordered_map<std::string, int> headerCount;
    std::regex headerPattern(R"(^\s*\[\s*([^\]]+?)\s*\]\s*$)");
    std::string line;
    int modified = 0;
    
    while (std::getline(infile, line)) {
        std::smatch matches;

        if (std::regex_match(line, matches, headerPattern)) {
            std::string headerName = matches[1].str();
            int count = ++headerCount[headerName];
            
            if (count > 1) {
                std::string newHeader = "[" + headerName + std::to_string(count) + "]";
                std::cout << "Modified: '" << line << "' -> '" << newHeader << "'" << std::endl;
                outfile << newHeader << std::endl;
                modified++;
            } else {
                outfile << line << std::endl;
            }
        } else {
            outfile << line << std::endl;
        }
    }
    
    infile.close();
    outfile.close();
    
    std::cout << "Processing completed! Modified " << modified << " headers" << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    
    return 0;
}