#pragma once
#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <sdsl/suffix_arrays.hpp>
#include "../include/k_factor_tree/k_factor_tree.h"

std::vector<std::pair<size_t, size_t>> construct_with_double_fm(const std::string& text, int tau_u, int lambda, const std::string& delimiter) {
    using fm_index_type = sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, 512, 1024> ;
    using fm_size_type = typename fm_index_type::size_type;

    fm_index_type forward_fm, backward_fm;
    sdsl::construct_im(forward_fm, std::string(text.rbegin(), text.rend()), 1);
    sdsl::construct_im(backward_fm, text, 1);

    std::vector<std::pair<size_t, size_t>> minimun_factors;
    std::unordered_set<char> delimiters;
    for (char c : delimiter) { delimiters.insert(c); }

    // boundary of the interval
    fm_size_type lb = 0, rb = forward_fm.size() - 1;
    for (size_t l = 0, r = 0; r < text.size(); r++) {
        if (delimiters.find(text[r]) != delimiters.end()) {
            l = r + 1;
            lb = 0;
            rb = forward_fm.size() - 1;
            continue;
        }
        
        auto cnt = sdsl::backward_search(forward_fm, lb, rb, text[r], lb, rb);

        if (cnt <= tau_u) {
            fm_size_type sub_lb = 0, sub_rb = backward_fm.size() - 1;
            if (r + 1 > l + lambda) { l = r - lambda + 1; } // l = max(l, r-lambda+1)
            for (size_t i = r; i >= l; --i){
                auto subCnt = sdsl::backward_search(backward_fm, sub_lb, sub_rb, text[i], sub_lb, sub_rb);
                if (subCnt <= tau_u) {
                    minimun_factors.emplace_back(i, r);
                    l = i + 1;
                    break;
                }
            }
            r = l - 1;
            lb = 0;
            rb = forward_fm.size() - 1;
        }
    }

    return minimun_factors;
}

TEST(gen_mf_test, compare_with_double_bwt) {
    size_t tau_l = 0;
    size_t tau_u = 6;
    size_t lambda = 20;
    //std::unordered_set<char> delimiters {'N'};
    std::string delimiter {"N"};

    std::string inputTextPath = "/Users/zhaozhewei/Desktop/exp/pseudo-real/dna/rawText"; // TODO: hard code
    //std::string inputTextPath = "/Users/zhaozhewei/Desktop/multi-layer-figiss/dataset/10MB_200People.txt"; // TODO: hard code
    std::ifstream inputFile(inputTextPath);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open input file " << inputTextPath << std::endl;
        return ;
    }

    std::stringstream buffer;
    buffer << inputFile.rdbuf();
    std::string intputContent = buffer.str();

    std::vector<std::pair<size_t, size_t>> mf1, mf2;
    std::cout << "k-factor-tree processing\n";
    k_factor_tree ksf(intputContent, lambda, tau_l, tau_u, delimiter);
    mf1 = ksf.get_min_factors();
    std::cout << "double BWT processing\n";
    mf2 = construct_with_double_fm(intputContent, tau_u, lambda, delimiter);
    ASSERT_EQ(mf1, mf2);
}