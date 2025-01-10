#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <sstream>
#include "../include/k_factor_tree/k_factor_tree.h"
#include "../include/tau_lambda_index.h"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] 
                  << " [input index path] [input pattern path] [output results path]" << std::endl;
        return 1;
    }

    std::string inputIndexPath = argv[1];
    std::string inputPatternPath = argv[2];
    std::string outputResultPath = argv[3];

    std::ifstream inputIndex(inputIndexPath);
    if (!inputIndex.is_open()) {
        std::cerr << "Error: Could not open input index file " << inputIndexPath << std::endl;
        return 1;
    }

    std::ifstream inputPatternFile(inputPatternPath);
    if (!inputPatternFile.is_open()) {
        std::cerr << "Error: Could not open input pattern file " << inputPatternPath << std::endl;
        return 1;
    }

    std::ofstream outputResultFile(outputResultPath);
    if (!outputResultFile.is_open()) {
        std::cerr << "Error: Could not open output result file " << outputResultPath << std::endl;
        return 1;
    }

    // create a tau-lambda-index instance and operate it
    //std::chrono::steady_clock::time_point t1, t2;
    tau_lambda_index* idx = new tau_lambda_index();
    idx->load(inputIndex);
    //t1 = std::chrono::steady_clock::now();
    idx->locate(inputPatternFile, outputResultFile);
    //t2 = std::chrono::steady_clock::now();

    inputIndex.close();
    inputPatternFile.close();
    outputResultFile.close();

    std::cout << "tau-lambda-index location queries completed successfully!" << std::endl;
    return 0;
}
