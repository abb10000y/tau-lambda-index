#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../include/k_factor_tree/k_factor_tree.h"

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] 
                  << " [input text path] [output mf path] [tau_l] [tau_u] [lambda] [delimiter/terminal symbols]" << std::endl;
        return 1;
    }

    std::string inputTextPath = argv[1];
    std::string outputMfPath = argv[2];
    uint64_t tau_l = std::stod(argv[3]);
    uint64_t tau_u = std::stod(argv[4]);
    uint64_t lambda = std::stod(argv[5]);
    std::string delimiter = argv[6];

    std::ifstream inputFile(inputTextPath);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open input file " << inputTextPath << std::endl;
        return 1;
    }

    std::ofstream outputFile(outputMfPath);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputMfPath << std::endl;
        return 1;
    }

    std::string intputContent;
    inputFile >> intputContent;

    k_factor_tree ksf(intputContent, lambda, tau_l, tau_u, delimiter);
    ksf.Serialize_min_factors(outputFile);

    inputFile.close();
    outputFile.close();

    std::cout << "Minimal factors identification completed successfully!" << std::endl;
    return 0;
}
