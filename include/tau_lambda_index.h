#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <exception>
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
#include "self_indexes/r-index/internal/r_index.hpp"
#include "self_indexes/r-index/internal/utils.hpp"
#include "self_indexes/LZ/src/static_selfindex.h"
#include "self_indexes/LZ/src/static_selfindex_lz77.h"
#include "self_indexes/LMS-based-self-index-LPG_grid/include/lpg/lpg_index.hpp"
#include "self_indexes/LMS-based-self-index-LPG_grid/third-party/CLI11.hpp"


class tau_lambda_index {
public:
    tau_lambda_index(){}
    tau_lambda_index(std::string text_path, std::string mf_path, size_t self_index_type);
    tau_lambda_index(std::string text_path, std::string mf_path, std::string index_path, size_t self_index_type);
    void serialize(std::ofstream &out);
    void load(std::ifstream &in, std::string inputIndexPath);
    void load_min_factors(std::ifstream &in);
    void locate(std::ifstream &in, std::ofstream &out);
    bool locate_xbwt(std::string &pattern);
    std::vector<ulint> locate_r_index(std::string &pattern);
    double get_masked_ratio() { return masked_ratio; }
    void log(std::ofstream& out) {
        out << "tau_l, tau_u, lambda: " << tau_l << ", " << tau_u << ", " << lambda << "\n";
        if (self_index_type == 1) {
            out << "bwt_runs: " << r_index->number_of_runs() << "\n";
        } else if (self_index_type == 2) {
            out << "number_of_phrases: " << "PLEASE CHECK THE CONSTRUCTION LOG" << "\n";
        }
        out << "masked_ratio: " << masked_ratio << "\n" ;
    }

private:
    size_t self_index_type, tau_l, tau_u, lambda;
    double masked_ratio;
    XBWT* xbwt {new XBWT()}; // TODO: memory leak
    SymbolTable symbol_table_;
    std::vector<std::pair<size_t, size_t>> min_factors; // TODO: change to local variable?
    //std::unordered_set<char> delimiters; // TODO: no need?
    uint64_t t_symbol = static_cast<uint64_t>(1); // terminate symbols

    // self_indexes
    ri::r_index<> *r_index {new ri::r_index<>()};
    lz77index::static_selfindex* lz77;
    lpg_index* lms;

    void build_XBWT(const std::string &text);
    void gen_masked_text(const std::string &text, std::string &masked_text);
    std::vector<uint64_t> _locate(std::string &pattern);
};

void tau_lambda_index::gen_masked_text(const std::string &text, std::string &masked_text) {
    std::vector<std::pair<size_t, size_t>> masked_notations;
    size_t n = text.size();
    // mark which parts should be kept
    if (min_factors.size() > 0) {
        size_t start = 0, end = n - 1;
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
    }

    // calculate the masked_ratio, (# of masked characters) / |text|
    size_t cnt = 0;
    for (auto v : masked_notations) {
        cnt += std::get<1>(v) - std::get<0>(v) + 1;
    }
    masked_ratio = 1.0 - 1.0 * cnt / n;

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

void tau_lambda_index::build_XBWT(const std::string &text) {
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
    in >> tau_l >> tau_u >> lambda >> tmp >> n;
    // for (auto c : tmp) { delimiters.insert(c); }
    while (n > 0) {
        size_t a, b;
        in >> a >> b;
        min_factors.push_back({a, b});
        n--;
    }
}

// Constructor for r-index and LMS
tau_lambda_index::tau_lambda_index(std::string text_path, std::string mf_path, size_t self_index_type): self_index_type(self_index_type) {
    std::ifstream text_in(text_path);
    if (!text_in.is_open()) {
        std::cerr << "Error: Could not open input text file " << text_path << std::endl;
        return;
    }
    std::stringstream buffer;
    buffer << text_in.rdbuf();
    std::string text = buffer.str();
    text_in.close();

    std::ifstream mf_in(mf_path);
    if (!mf_in.is_open()) {
        std::cerr << "Error: Could not open input minimal_factors file " << mf_path << std::endl;
        return;
    }
    load_min_factors(mf_in);
    mf_in.close();

    build_XBWT(text);
    std::string maskedTextPath = "tmpMaskedText";
    std::string maskedText;
    std::ofstream tmpMaskedTextfile(maskedTextPath);
    if (!tmpMaskedTextfile.is_open()) {
        std::cerr << "Error: Could not open output tmpMaskedText file " << maskedTextPath << std::endl;
        return;
    }

    gen_masked_text(text, maskedText);
    tmpMaskedTextfile << maskedText;
    tmpMaskedTextfile.close();

    // TODO: generate the self-index
    if (self_index_type == 1) {
        std::string input;
        std::ifstream fs(maskedTextPath);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
        r_index = new ri::r_index<>(input, true);
    } else if (self_index_type == 3) {
        std::string LMSTemp = "LMS_temp";
        lms = new lpg_index(text_path, LMSTemp, 1, 0.5);
    }

    if (std::remove(maskedTextPath.c_str()) != 0) {
        std::cerr << "Not able to delete the masked text tmp file" << std::endl;
        return;
    }
}

// Constructor for LZ77
tau_lambda_index::tau_lambda_index(std::string text_path, std::string mf_path, std::string index_path, size_t self_index_type): self_index_type(self_index_type) {
    std::ifstream text_in(text_path);
    if (!text_in.is_open()) {
        std::cerr << "Error: Could not open input text file " << text_path << std::endl;
        return;
    }
    std::stringstream buffer;
    buffer << text_in.rdbuf();
    std::string text = buffer.str();
    text_in.close();

    std::ifstream mf_in(mf_path);
    if (!mf_in.is_open()) {
        std::cerr << "Error: Could not open input minimal_factors file " << mf_path << std::endl;
        return;
    }
    load_min_factors(mf_in);
    mf_in.close();

    build_XBWT(text);
    std::string maskedTextPath = "tmpMaskedText";
    std::string maskedText;
    std::ofstream tmpMaskedTextfile(maskedTextPath);
    if (!tmpMaskedTextfile.is_open()) {
        std::cerr << "Error: Could not open output tmpMaskedText file " << maskedTextPath << std::endl;
        return;
    }

    gen_masked_text(text, maskedText);
    tmpMaskedTextfile << maskedText;
    tmpMaskedTextfile.close();

    if (self_index_type == 2) {
        unsigned char br=0;
        unsigned char bs=0;
        unsigned char ss=0;
        char *in = new char[maskedTextPath.size() + 1], *out = new char[index_path.size() + 1]; // 分配記憶體（+1 用於結尾的 '\0'）
        std::strcpy(in, maskedTextPath.c_str());
        std::strcpy(out, index_path.c_str());
        lz77 = lz77index::static_selfindex_lz77::build(in, out, br, bs, ss);
    }

    if (std::remove(maskedTextPath.c_str()) != 0) {
        std::cerr << "Not able to delete the masked text tmp file" << std::endl;
        return;
    }
}

void tau_lambda_index::serialize(std::ofstream &out) {
    sdsl::write_member(self_index_type, out);
    sdsl::write_member(tau_l, out);
    sdsl::write_member(tau_u, out);
    sdsl::write_member(lambda, out);
    sdsl::write_member(masked_ratio, out);
    symbol_table_.Serialize(out);
    xbwt->Serialize(out);
    if (self_index_type == 1) {
        r_index->serialize(out);
    } else if (self_index_type == 3) {
        lms->serialize(out, NULL, "");
    }
}

void tau_lambda_index::load(std::ifstream &in, std::string inputIndexPath) {
    sdsl::read_member(self_index_type, in);
    sdsl::read_member(tau_l, in);
    sdsl::read_member(tau_u, in);
    sdsl::read_member(lambda, in);
    sdsl::read_member(masked_ratio, in);
    symbol_table_.Load(in);
    xbwt->Load(in);
    if (self_index_type == 1) {
        r_index->load(in);
    } else if (self_index_type == 2) {
        std::string lz77Path = inputIndexPath + "_lz77";
        char *in = new char[lz77Path.size() + 1]; // 分配記憶體（+1 用於結尾的 '\0'）
        std::strcpy(in, lz77Path.c_str());
        lz77 = lz77index::static_selfindex::load(in);
    }
}

std::vector<uint64_t> tau_lambda_index::_locate(std::string &pattern) {
    vector<uint64_t> result;
    sdsl::int_vector<> pattern_int;
    pattern_int.width(64);
    pattern_int.resize(pattern.size());
    for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[pattern[i] + t_symbol]; }
    // auto [start_pattern, length] = xbwt->match_pos_in_pattern(pattern_int);
    
    if (pattern.length() <= lambda && xbwt->match_if_exist(pattern_int)) {
        if (self_index_type == 1) {
            result = r_index->locate_all_tau(pattern, tau_l);
        } else if (self_index_type == 2) {
            unsigned char *p = new unsigned char[pattern.size() + 1]; // 分配記憶體（+1 用於結尾的 '\0'）
            std::memcpy(p, pattern.c_str(), pattern.size() + 1);
            unsigned int nooc;
            std::vector<unsigned int> *tmp = lz77->locate(p, pattern.size(), &nooc);
            result.resize(tmp->size());
            for (size_t i = 0; i < tmp->size(); i++) { result[i] = static_cast<uint64_t>((*tmp)[i]); }
        }
    }

    return result;
}
void tau_lambda_index::locate(std::ifstream &in, std::ofstream &out) {
    std::chrono::steady_clock::time_point t1, t2;

    string header;
	std::getline(in, header);

	size_t n = get_number_of_patterns(header);
	size_t m = get_patterns_length(header);

    for (size_t i = 0; i < n; i++) {
        string pattern = string();

		for(size_t j = 0; j < m; j++){
			char c;
			in.get(c);
			pattern += c;
		}

        t1 = std::chrono::steady_clock::now();
        std::vector<uint64_t> results = _locate(pattern);
        t2 = std::chrono::steady_clock::now();
        std::sort(results.begin(), results.end());
        out << "pattern [" << (i+1) << "] -- " << results.size() <<  "occurrences\n";
        for (auto r : results) { out << r << "\t"; }
        out << "\n";
        out << "time consuming (us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "\n";
    }
}

// TO be deleted
bool tau_lambda_index::locate_xbwt(std::string &pattern) {
    sdsl::int_vector<> pattern_int;
    pattern_int.width(64);
    pattern_int.resize(pattern.size());
    for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[pattern[i] + t_symbol]; }
    // auto [start_pattern, length] = xbwt->match_pos_in_pattern(pattern_int);
        
    return xbwt->match_if_exist(pattern_int);
}

// TO be deleted
std::vector<ulint> tau_lambda_index::locate_r_index(std::string &pattern) {
    // vector<ulint> result;
    // result = r_index->locate_all_tau(pattern, tau_l, tau_u);
    // return result;
}