// an simple baseline index, built by fm-index in sdsl-lite

#pragma once
#include <sdsl/suffix_arrays.hpp>
#include <exception>

class Baseline_index {
public:
    Baseline_index(const std::string& text, int tao, int lambda = -1, bool activate_lambda = false);
    std::vector<size_t> search(const std::string& pattern, bool verbose = true);

private:
    const std::string& text;
    int tao;
    int lambda;
    bool activate_lambda;
};

Baseline_index::Baseline_index(const std::string& text, int tao, int lambda, bool activate_lambda)
    : text(text), tao(tao), lambda(lambda), activate_lambda(activate_lambda) {
    if (activate_lambda && lambda < 0) {
        std::cout << "Error: lambda is not set" << std::endl;
        throw std::runtime_error("lambda is not set"); // TODO more specific exception?
    }
}

std::vector<size_t> Baseline_index::search(const std::string& pattern, bool verbose) {
    if (activate_lambda && pattern.size() > lambda) {
        if (verbose) {
            std::cout << "Warning: the length of the pattern is larger than lambda" << std::endl;
        }
        return std::vector<size_t>();
    }

    sdsl::csa_sada<> csa;
    sdsl::construct_im(csa, text, 1);
    auto csa_result = sdsl::locate(csa, pattern.begin(), pattern.end());

    if (csa_result.size() > tao) {
        return std::vector<size_t>();
    }

    if (csa_result.size() == 0) {
        return std::vector<size_t>();
    }
    
    return std::vector<size_t>(csa_result.begin(), csa_result.end());
}