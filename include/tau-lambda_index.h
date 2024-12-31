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
#include "../r-index/internal/r_index.hpp"


class tau-lambda_index {
public:
    tau-lambda_index(){}
    tau-lambda_index(std::string text_path, int tau_lower, int tau_upper, int lambda);
    void serialize(std::ofstream &out_path);
    void load(std::ifstream &in_path);
    std::vector<ulint> locate(std::string &pattern);
    bool locate_xbwt(std::string &pattern);
    std::vector<ulint> locate_r_index(std::string &pattern);
    size_t get_bwt_runs() { return r_index->number_of_runs(); }

private:
    bool is_raw_r_index = false; // TODO: judgement not finished
    XBWT* xbwt {new XBWT()};
    SymbolTable symbol_table_;
    ri::r_index<> *r_index {new ri::r_index<>()};
    uint64_t t_symbol = static_cast<uint64_t>(1); // terminate symbols
    size_t tau_lower, tau_upper, lambda;
};

tau-lambda_index::tau-lambda_index(std::string text_path, int tau_lower, int tau_upper, int lambda)
    : tau_lower(tau_lower), tau_upper(tau_upper), lambda(lambda)
{
    ri::r_index<> *r_index_raw, *r_index_masked;

    if (tau_lower > tau_upper || tau_lower < 0 || lambda < 0) {
        throw std::runtime_error("error input parameter"); // TODO more specific exception?
    }
    if (!std::filesystem::exists(text_path)) {
        throw std::runtime_error("input file doesn't exist"); // TODO more specific exception?
    }

    std::string text;
    std::ifstream text_in(text_path);
    std::getline(text_in, text);
    text_in.close();

    // gen r_index_raw
    {
        std::string input;
        std::ifstream fs(text_path);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
        r_index_raw = new ri::r_index<>(input, true);
    }

    // find min_factor from k_factor_tree
    std::vector<std::pair<size_t, size_t>> mf;
    std::string masked_text;
    {
        k_factor_tree ksf (text, lambda, tau_lower, tau_upper);
        mf = ksf.get_min_factors();
        ksf.gen_masked_text(text, masked_text, lambda);
    }

    // build AC_automaon for min_factors
    std::vector<uint64_t> concated_mf; // each substring is separated by terminal_symbol (t_symbol)
    for (auto [b, e] : mf) {
        for (size_t i = e; i >= b; i--) { concated_mf.push_back(text[i] + t_symbol); } // reverse insert
        concated_mf.push_back(t_symbol);
    }
    sdsl::int_vector<> concated_mf_int_vector;
    // size_t usedBit = 9; // TODO: use variable
    concated_mf_int_vector.width(64);
    concated_mf_int_vector.resize(concated_mf.size());
    for (size_t i = 0; i < concated_mf.size(); i++) { concated_mf_int_vector[i] = concated_mf[i]; }
    sdsl::append_zero_symbol(concated_mf_int_vector);
    symbol_table_ = decltype(symbol_table_)(concated_mf_int_vector);
    for (size_t i = 0; i < concated_mf_int_vector.size(); i++) { concated_mf_int_vector[i] = symbol_table_[concated_mf_int_vector[i]]; }
    xbwt->insert(concated_mf_int_vector, symbol_table_.GetEffectiveAlphabetWidth(), symbol_table_.GetAlphabetSize());

    // gen r_index_masked
    {
        std::string tmp_path = "/mnt/f/alg/git/multi-layer-figiss/experiments/masked_tmp"; // TODO: hard code
        std::ofstream out(tmp_path);
        out << masked_text;
        out.close();

        std::string input;
        std::ifstream fs(tmp_path);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
        r_index_masked = new ri::r_index<>(input, true);
        std::remove(tmp_path.c_str());
    }

    // serialize and reload all (test only)
    /*
    {
        std::string path_1 = "/mnt/f/alg/git/multi-layer-figiss/experiments/masked_r_index"; // TODO: hard code
        std::string path_2 = "/mnt/f/alg/git/multi-layer-figiss/experiments/masked_xbwt"; // TODO: hard code
        std::string path_3 = "/mnt/f/alg/git/multi-layer-figiss/experiments/raw_r_index"; // TODO: hard code

        std::ofstream out_1(path_1);
        std::ofstream out_2(path_2);
        std::ofstream out_3(path_3);

        r_index_masked.serialize(out_1);
        xbwt->Serialize(out_2);
        symbol_table_.Serialize(out_2);
        r_index_raw.serialize(out_3);

        out_1.close();
        out_2.close();
        out_3.close();

        XBWT* xbwt_2 {new XBWT()};
        SymbolTable symbol_table_2;
        ri::r_index<> r_index_raw_2, r_index_masked_2;

        std::ifstream in_1(path_1);
        std::ifstream in_2(path_2);
        std::ifstream in_3(path_3);

        r_index_masked_2.load(in_1);
        xbwt_2->Load(in_2);
        symbol_table_2.Load(in_2);
        r_index_raw_2.load(in_3);

        in_1.close();
        in_2.close();
        in_3.close();

        cout << "size of masked r-index: " << (sizeof(r_index_masked_2) + sizeof(*xbwt) + sizeof(symbol_table_2)) << " (byte)\n";
        cout << "size of masked r-index (r-index): " << sizeof(r_index_masked_2) << " (byte)\n";
        cout << "size of masked r-index (xbwt): " << (sizeof(*xbwt) + sizeof(symbol_table_2)) << " (byte)\n";
        cout << "size of raw r-index: " << sizeof(r_index_raw_2) << " (byte)\n";
    }
    */

    r_index = r_index_masked;
}

void tau-lambda_index::serialize(std::ofstream &out_path) {
    sdsl::write_member(is_raw_r_index, out_path);
    sdsl::write_member(tau_lower, out_path);
    sdsl::write_member(tau_upper, out_path);
    sdsl::write_member(lambda, out_path);
    if (is_raw_r_index) {
        r_index->serialize(out_path);
    } else {
        symbol_table_.Serialize(out_path);
        xbwt->Serialize(out_path);
        r_index->serialize(out_path);
    }
}

void tau-lambda_index::load(std::ifstream &in_path) {
    sdsl::read_member(is_raw_r_index, in_path);
    sdsl::read_member(tau_lower, in_path);
    sdsl::read_member(tau_upper, in_path);
    sdsl::read_member(lambda, in_path);
    if (is_raw_r_index) {
        r_index->load(in_path);
    } else {
        symbol_table_.Load(in_path);
        xbwt->Load(in_path);
        r_index->load(in_path);
    }
}

std::vector<ulint> tau-lambda_index::locate(std::string &pattern) {
    vector<ulint> result;
    if (is_raw_r_index) {
        result = r_index->locate_all_tau(pattern, tau_lower, tau_upper);
    } else {
        sdsl::int_vector<> pattern_int;
        pattern_int.width(64);
        pattern_int.resize(pattern.size());
        for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[pattern[i] + t_symbol]; }
        // auto [start_pattern, length] = xbwt->match_pos_in_pattern(pattern_int);
        
        if (xbwt->match_if_exist(pattern_int)) {
            result = r_index->locate_all_tau(pattern, tau_lower, tau_upper);
        }
    }
    return result;
}

bool tau-lambda_index::locate_xbwt(std::string &pattern) {
    sdsl::int_vector<> pattern_int;
    pattern_int.width(64);
    pattern_int.resize(pattern.size());
    for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[pattern[i] + t_symbol]; }
    // auto [start_pattern, length] = xbwt->match_pos_in_pattern(pattern_int);
        
    return xbwt->match_if_exist(pattern_int);
}

std::vector<ulint> tau-lambda_index::locate_r_index(std::string &pattern) {
    vector<ulint> result;
    result = r_index->locate_all_tau(pattern, tau_lower, tau_upper);
    return result;
}