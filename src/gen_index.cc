#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "../include/k_factor_tree/k_factor_tree.h"
#include "../include/tau_lambda_index.h"

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] 
                  << " [input text path] [mf path] [output index path] [self-index type]" << std::endl;
        return 1;
    }

    std::string inputTextPath = argv[1];
    std::string inputMfPath = argv[2];
    std::string outputIndexPath = argv[3];
    uint64_t selfIndexType = std::stod(argv[4]);

    std::ofstream outputFile(outputIndexPath); // TODO: move to the tau-lambda-index header file
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputIndexPath << std::endl;
        return 1;
    }

    // create a tau-lambda-index instance and operate it

    outputFile.close();

    std::cout << "tau-lambda-index completed successfully!" << std::endl;
    return 0;
}
