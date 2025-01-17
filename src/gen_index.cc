#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <sstream>
#include <chrono>
#include "../include/k_factor_tree/k_factor_tree.h"
#include "../include/tau_lambda_index.h"

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] 
                  << " [input text path] [mf path] [output index path] [self-index type] [log path](optional)" << std::endl;
        return 1;
    }

    std::string inputTextPath = argv[1];
    std::string inputMfPath = argv[2];
    std::string outputIndexPath = argv[3];
    // uint64_t selfIndexType = std::stod(argv[4]);
    index_types index_type = static_cast<index_types>(std::stod(argv[4]));
    std::string outputLogPath;
    if (argc == 6) { outputLogPath = argv[5]; }

    std::ofstream outputFile(outputIndexPath);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputIndexPath << std::endl;
        return 1;
    }

    std::ofstream outputLogFile(outputLogPath);
    if (argc == 6 && !outputLogFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputLogPath << std::endl;
        return 1;
    }

    // create a tau-lambda-index instance and operate it
    tau_lambda_index* idx;
    std::chrono::steady_clock::time_point t1, t2;
    t1 = std::chrono::steady_clock::now();
    if (index_type == index_types::r_index_type || index_type == index_types::LMS_type || index_type == index_types::old_tau_lambda_type) {
        idx = new tau_lambda_index(inputTextPath, inputMfPath, index_type);
    } else if (index_type == index_types::lz77_type) {
        std::string lz77OutPath= outputIndexPath + "_lz77";
        idx = new tau_lambda_index(inputTextPath, inputMfPath, lz77OutPath, index_type);
    } else {
        throw std::runtime_error("Invalid selfIndexType\n");
    }
    idx->serialize(outputFile);
    t2 = std::chrono::steady_clock::now();

    if (argc == 6) {
        std::filesystem::path path(inputTextPath);
        std::string fileName = path.stem().string();
        outputLogFile << "[" << fileName << "]\n";
        idx->log(outputLogFile);
        outputLogFile << "consumed time (us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "\n";
        outputLogFile.close();
    }

    outputFile.close();

    std::cout << "tau-lambda-index completed successfully!" << std::endl;
    return 0;
}
