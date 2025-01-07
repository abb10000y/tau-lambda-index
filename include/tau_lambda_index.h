#pragma once

#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <string_view>
#include <exception>
#include <filesystem>
#include <unordered_set>
#include <chrono>
#include <sys/resource.h>
#include <sys/time.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include "xbwt/xbwt.h"
#include "symbol_table/symbol_table.h"
#include "k_factor_tree/k_factor_tree.h"
#include "util/utility.h"


class tau_lambda_index {
public:
    tau_lambda_index(){}
    tau_lambda_index(std::string text_path, std::string mf_path);
    void serialize(std::ofstream &out_path);
    void load(std::ifstream &in_path);
    void load_min_factors(std::ifstream &in);
    std::vector<ulint> locate(std::string &pattern);
    bool locate_xbwt(std::string &pattern);
    std::vector<ulint> locate_r_index(std::string &pattern);
    double get_coverage_rate() { return coverage_rate; }

private:
    // bool is_raw_r_index = false; // TODO: judgement not finished
    std::string text;
    XBWT* xbwt {new XBWT()}; // TODO: memory leak
    SymbolTable symbol_table_;
    uint64_t t_symbol = static_cast<uint64_t>(1); // terminate symbols
    size_t tau_l, tau_u, lambda;
    std::vector<std::pair<size_t, size_t>> min_factors;
    std::unordered_set<char> delimiters;
    double coverage_rate;

    void build_XBWT();
    void gen_masked_text(const std::string &text, std::string &masked_text);
};

void tau_lambda_index::gen_masked_text(const std::string &text, std::string &masked_text) {
    if (min_factors.size() == 0) { throw std::invalid_argument("min_factors is empty"); }
    std::vector<std::pair<size_t, size_t>> masked_notations;

    // mark which parts should be kept
    size_t n = text.size(), start = 0, end = n - 1;
    if (std::get<1>(min_factors[0]) + 1 > lambda) { start = std::get<1>(min_factors[0]) + 1 - lambda; }
    if (std::get<0>(min_factors[0]) + lambda - 1 < n) { end = std::get<0>(min_factors[0]) + lambda - 1; }
    for (size_t i = 1, m = min_factors.size(); i < m; i++) {
        size_t next_start = 0, next_end = n - 1;
        if (std::get<1>(min_factors[i]) + 1 > lambda) { next_start = std::get<1>(min_factors[i]) + 1 - lambda; }
        if (std::get<0>(min_factors[i]) + lambda - 1 < n) { next_end = std::get<0>(min_factors[i]) + lambda - 1; }
        if (next_start <= end) { end = next_end; }
        else {
            masked_notations.push_back({start, end});
            start = next_start;
            end = next_end;
        }
    }
    masked_notations.push_back({start, end});

    // calculate the coverage_rate, (# of masked characters) / |text|
    size_t cnt = 0;
    for (auto v : masked_notations) {
        cnt += std::get<1>(v) - std::get<0>(v) + 1;
    }
    coverage_rate = 1.0 * cnt / n;

    // generate the masked text
    size_t masked_symbol = 255; // TODO: hard code
    if (lambda == 0) { lambda = text.size(); }
    if (masked_notations.size() == 0) { throw std::invalid_argument("masked_notations is empty"); } // TODO: if this necessary?

    masked_text.assign(n, masked_symbol);
    for (auto v : masked_notations) {
        for (size_t i = std::get<0>(v); i <= std::get<1>(v); i++) {
            masked_text[i] = text[i];
        }
    }
}

void tau_lambda_index::build_XBWT() {
    std::vector<uint64_t> concated_mf; // each substring is separated by terminal_symbol (t_symbol)
    for (auto [b, e] : min_factors) {
        for (size_t i = e; i >= b; i--) { concated_mf.push_back(text[i] + t_symbol); } // reverse insert
        concated_mf.push_back(t_symbol);
    }
    sdsl::int_vector<> concated_mf_int_vector;
    concated_mf_int_vector.width(64); // TODO: change to the smallest bits usage
    concated_mf_int_vector.resize(concated_mf.size());
    for (size_t i = 0; i < concated_mf.size(); i++) { concated_mf_int_vector[i] = concated_mf[i]; }
    sdsl::append_zero_symbol(concated_mf_int_vector);
    symbol_table_ = decltype(symbol_table_)(concated_mf_int_vector);
    for (size_t i = 0; i < concated_mf_int_vector.size(); i++) { concated_mf_int_vector[i] = symbol_table_[concated_mf_int_vector[i]]; }
    xbwt->insert(concated_mf_int_vector, symbol_table_.GetEffectiveAlphabetWidth(), symbol_table_.GetAlphabetSize());
}

void tau_lambda_index::load_min_factors(std::ifstream &in) {
    std::string tmp;
    size_t n;
    in >> lambda >> tau_l >> tau_u >> tmp >> n;
    for (auto c : tmp) { delimiters.insert(c); }
    while (n > 0) {
        size_t a, b;
        in >> a >> b;
        min_factors.push_back({a, b});
        n--;
    }
}

tau_lambda_index::tau_lambda_index(std::string text_path, std::string mf_path) {
    std::ifstream text_in(text_path);
    if (!text_in.is_open()) {
        std::cerr << "Error: Could not open input text file " << text_path << std::endl;
        return;
    }
    text_in >> text;
    text_in.close();

    std::ifstream mf_in(mf_path);
    if (!mf_in.is_open()) {
        std::cerr << "Error: Could not open input minimal_factors file " << mf_path << std::endl;
        return;
    }
    load_min_factors(mf_in);
    mf_in.close();

    build_XBWT();
    std::string maskedTextPath = "tmpMaskedText";
    std::string maskedText;
    std::ofstream tmpMaskedTextfile(maskedTextPath);
    if (!tmpMaskedTextfile.is_open()) {
        std::cerr << "Error: Could not open output tmpMaskedText file " << inputMfPath << std::endl;
        return;
    }

    gen_masked_text(text, maskedText);
    tmpMaskedTextfile << maskedText;
    tmpMaskedTextfile.close();

    // TODO: generate the self-index

    if (std::remove(maskedTextPath.c_str()) != 0) {
        std::cerr << "Not able to delete the masked text tmp file" << std::endl;
        return 1;
    }
}

void tau_lambda_index::serialize(std::ofstream &out_path) {
    // sdsl::write_member(is_raw_r_index, out_path);
    // sdsl::write_member(tau_l, out_path);
    // sdsl::write_member(tau_u, out_path);
    // sdsl::write_member(lambda, out_path);
    // if (is_raw_r_index) {
    //     r_index->serialize(out_path);
    // } else {
    //     symbol_table_.Serialize(out_path);
    //     xbwt->Serialize(out_path);
    //     r_index->serialize(out_path);
    // }
}

void tau_lambda_index::load(std::ifstream &in_path) {
    // sdsl::read_member(is_raw_r_index, in_path);
    // sdsl::read_member(tau_l, in_path);
    // sdsl::read_member(tau_u, in_path);
    // sdsl::read_member(lambda, in_path);
    // if (is_raw_r_index) {
    //     r_index->load(in_path);
    // } else {
    //     symbol_table_.Load(in_path);
    //     xbwt->Load(in_path);
    //     r_index->load(in_path);
    // }
}

std::vector<ulint> tau_lambda_index::locate(std::string &pattern) {
    // vector<ulint> result;
    // if (is_raw_r_index) {
    //     result = r_index->locate_all_tau(pattern, tau_l, tau_u);
    // } else {
    //     sdsl::int_vector<> pattern_int;
    //     pattern_int.width(64);
    //     pattern_int.resize(pattern.size());
    //     for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[pattern[i] + t_symbol]; }
    //     // auto [start_pattern, length] = xbwt->match_pos_in_pattern(pattern_int);
        
    //     if (xbwt->match_if_exist(pattern_int)) {
    //         result = r_index->locate_all_tau(pattern, tau_l, tau_u);
    //     }
    // }
    // return result;
}

bool tau_lambda_index::locate_xbwt(std::string &pattern) {
    sdsl::int_vector<> pattern_int;
    pattern_int.width(64);
    pattern_int.resize(pattern.size());
    for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[pattern[i] + t_symbol]; }
    // auto [start_pattern, length] = xbwt->match_pos_in_pattern(pattern_int);
        
    return xbwt->match_if_exist(pattern_int);
}

std::vector<ulint> tau_lambda_index::locate_r_index(std::string &pattern) {
    // vector<ulint> result;
    // result = r_index->locate_all_tau(pattern, tau_l, tau_u);
    // return result;
}