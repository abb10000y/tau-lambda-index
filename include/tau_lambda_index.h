#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
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
#include "self_indexes/multi-layer-figiss/src/components/Index.h"
#include "self_indexes/hybrid/hybrid.h"
#include "compact_suffix_trie.h"

#ifndef ulong
#define ulong unsigned long
#endif

enum class index_types {
    r_index_type = 1,
    lz77_type = 2,
    LMS_type = 3,
    old_tau_lambda_type = 4,
    compact_suffix_trie = 5,
    hybrid = 6
};

class tau_lambda_index {
public:
    tau_lambda_index(){}
    tau_lambda_index(std::string &mf_path);
    tau_lambda_index(std::string &mf_path, index_types index_type);
    tau_lambda_index(std::string &mf_path, std::string &index_path, index_types index_type);
    tau_lambda_index(std::string &mf_path, std::string &index_path, std::string &outputFolder, index_types index_type, uint M);
    void serialize(std::ofstream &out);
    void load(std::ifstream &in, std::string inputIndexPath, bool xbwt_only = false);
    void locate(std::ifstream &in, std::ofstream &out, bool xbwt_only = false);
    double get_masked_ratio() { return masked_ratio; }
    void log(std::ofstream& out) {
        out << "[" << inputTextPath << "]\n";
        out << "tau_l, tau_u, lambda: " << tau_l << ", " << tau_u << ", " << lambda << "\n";
        if (index_type == index_types::r_index_type) {
            out << "bwt_runs: " << r_index->number_of_runs() << "\n";
        } else if (index_type == index_types::lz77_type) {
            out << "number_of_phrases: " << lz77->z << "\n";
        } else if (index_type == index_types::LMS_type) {
            out << "grammar_size: " << lms.grammar_tree.get_grammar_size() << "\n";
        }
        out << "masked_ratio: " << masked_ratio << "\n" ;
    }
    void gen_patterns(std::string &path, size_t cnt) {
        std::ofstream out_with(path + "_with"), out_without(path + "_without");
        out_with << "# number=" << cnt << " length=" << lambda << " file=" << inputTextPath
                 << " tau_l=" << tau_l << " tau_u=" << tau_u << " lambda=" << lambda << " \n";
        out_without << "# number=" << cnt << " length=" << lambda << " file=" << inputTextPath
                 << " tau_l=" << tau_l << " tau_u=" << tau_u << " lambda=" << lambda << " \n";

        std::vector<std::pair<size_t, size_t>> with_candidates, without_candidates;

        size_t n = text.size();
        // start point candidates for query patterns with min_factor
        {
            size_t s = min_factors[0].first, e = s;
            // j-lambda+1>=0
            if (min_factors[0].second + 1 >= lambda) { s = min_factors[0].second + 1 - lambda; }
            for (auto [i, j] : min_factors) {
                // j-lambda+1<=e+1
                if (e + lambda >= j) { e = i; }
                else {
                    with_candidates.push_back({s, e});
                    s = j + 1 - lambda;
                    e = i;
                }
            }
            with_candidates.push_back({s, e});
        }
        // start point candidates for query patterns without min_factor
        {
            size_t s = 0, idx = 0;
            // e-lambda+1>=0
            if (with_candidates[idx].first == 0) {
                s = with_candidates[idx].second + 1;
                idx++;
            }
            for (size_t end = with_candidates.size(); idx < end; idx++) {
                auto [i, j] = with_candidates[idx];
                // i-1 <= n-1
                without_candidates.push_back({s, i - 1});
                s = j + 1;
            }
            if (s < n) { without_candidates.push_back({s, n - 1}); }
        }
        
        // random setting
        std::random_device rd;
        std::mt19937 gen(rd());
        size_t n_with = with_candidates.size(), n_without = without_candidates.size();

        // with min_factors
        for (size_t i = 0; i < cnt; i++) {
            std::uniform_int_distribution<> dis_a(0, n_with - 1);
            size_t a = dis_a(gen), len = with_candidates[a].second - with_candidates[a].first + 1;
            std::uniform_int_distribution<> dis_b(0, len - 1);
            size_t b = dis_b(gen);
            out_with << text.substr(with_candidates[a].first + b, lambda);
        }

        // without min_factors
        for (size_t i = 0; i < cnt; i++) {
            std::uniform_int_distribution<> dis_a(0, n_without - 1);
            size_t a = dis_a(gen), len = without_candidates[a].second - without_candidates[a].first + 1;
            std::uniform_int_distribution<> dis_b(0, len - 1);
            size_t b = dis_b(gen);
            out_without << text.substr(without_candidates[a].first + b, lambda);
        }

        out_with.close();
        out_without.close();
    }

private:
    size_t tau_l, tau_u, lambda;
    char delimiter = '#';
    index_types index_type;
    std::string inputTextPath, text;
    double masked_ratio;
    std::unique_ptr<XBWT> xbwt = std::make_unique<XBWT>();
    SymbolTable symbol_table_;
    std::vector<std::pair<size_t, size_t>> min_factors; // TODO: change to local variable?
    //std::unordered_set<char> delimiters; // TODO: no need?
    uint64_t t_symbol = static_cast<uint64_t>(1); // terminate symbols
    bool is_original_index = false;

    // underlying_indexes
    ri::r_index<> *r_index; 
    lz77index::static_selfindex_lz77* lz77;
    lpg_index lms;
    Index* old_tau_lambda;
    compact_suffix_trie* location_trie;
    hybrid* hybrid_index;

    void load_min_factors(std::string &mf_path);
    void build_XBWT(const std::string &text);
    void gen_masked_text(const std::string &text, std::string &masked_text);
    void _locate(std::string &pattern, std::vector<uint64_t> &results, bool xbwt_only = false);
    void _locate_original_index(std::string &pattern, std::vector<uint64_t> &results, bool xbwt_only = false);
    void load_Text(std::string &text, std::string &textPath) {
        std::ifstream text_in(textPath);
        if (!text_in.is_open()) {
            std::cerr << "Error: Could not open input text file " << textPath << std::endl;
            return;
        }
        std::stringstream buffer;
        buffer << text_in.rdbuf();
        text = buffer.str();
        text_in.close();
    }
    void output_Text(std::string &text, std::string &textPath) {
        std::ofstream out(textPath);
        if (!out.is_open()) {
            std::cerr << "Error: Could not open output file " << textPath << std::endl;
            return;
        }
        out << text;
        out.close();
    }

    // for generating query patterns
    void gen_covered_notations(std::vector<std::pair<size_t, size_t>> &covered_notations) {
        size_t n = text.size(), l_delimiter = 0, r_delimiter = 0;

        auto find_next_delimiter [&](size_t &l, size_t &r) -> void{
            if (r < n) {
                r = l;
                while (r < n && text[r] != delimiter) { r++; }
            }
        };

        // mark which parts should be kept
        if (min_factors.size() > 0) {
            find_next_delimiter(l_delimiter, r_delimiter);
            size_t start = 0, end = 0, 
            if (std::get<1>(min_factors[0]) + 2 > lambda) { start = std::get<1>(min_factors[0]) + 1 - lambda; } // TBD: some bugs here, depends on if T[0] is delimiter or not
            end = std::min(std::get<0>(min_factors[0]) + lambda - 1, r_delimiter - 1);
            for (size_t i = 1, m = min_factors.size(); i < m; i++) {
                size_t next_start = 0, next_end = r_delimiter - 1;
                if (std::get<0>(min_factors[i]) > r_delimiter) { find_next_delimiter(l_delimiter, r_delimiter); }
                if (std::get<1>(min_factors[i]) + 1 > lambda + l_delimiter) { next_start = std::get<1>(min_factors[i]) + 1 - lambda + l_delimiter; }
                else { next_start = l_delimiter + 1; } // TBD: some bugs here, depends on if T[0] is delimiter or not
                next_end = std::min(std::get<0>(min_factors[i]) + lambda - 1, r_delimiter - 1);
                if (next_start <= end) { end = next_end; }
                else {
                    covered_notations.push_back({start, end});
                    start = next_start;
                    end = next_end;
                }
            }
            covered_notations.push_back({start, end});
        }
    }
    size_t get_pattern_info(std::string keyword, std::string &header) {      
        ulint start_pos = header.find(keyword);
        if (start_pos == std::string::npos or start_pos+keyword.size()>=header.size())
            throw std::invalid_argument("Invalid query pattern format, can't find " + keyword + ". Please use gen_pattern to generate query patterns.");

        start_pos += keyword.size();

        ulint end_pos = header.substr(start_pos).find(" ");
        if (end_pos == std::string::npos)
            throw std::invalid_argument("Invalid query pattern format, can't find " + keyword + ". Please use gen_pattern to generate query patterns.");

        ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

        return n;
    }
};

void tau_lambda_index::gen_masked_text(const std::string &text, std::string &masked_text) {
    if (is_original_index) { masked_text = text; }
    else {
        std::vector<std::pair<size_t, size_t>> covered_notations;
        gen_covered_notations(covered_notations);

        // calculate the masked_ratio, (# of masked characters) / |text|
        size_t cnt = 0, n = text.size();
        for (auto v : covered_notations) {
            cnt += std::get<1>(v) - std::get<0>(v) + 1;
        }
        masked_ratio = 1.0 - 1.0 * cnt / n;

        // generate the masked text
        size_t masked_symbol = 255; // TODO: hard code
        masked_text.assign(n, masked_symbol);
        for (auto v : covered_notations) {
            for (size_t i = std::get<0>(v); i <= std::get<1>(v); i++) {
                masked_text[i] = text[i];
            }
        }
    }
}

void tau_lambda_index::build_XBWT(const std::string &text) {
    if (tau_u > 0) {
        std::vector<uint64_t> concated_mf; // each substring is separated by terminal_symbol (t_symbol)
        for (auto [b, e] : min_factors) {
            for (size_t i = e; i > b; i--) { concated_mf.push_back(static_cast<unsigned char>(text[i]) + t_symbol); } // reverse insert
            concated_mf.push_back(static_cast<unsigned char>(text[b]) + t_symbol); // prevent underflow for b == 0
            concated_mf.push_back(t_symbol);
        }
        sdsl::int_vector<> concated_mf_int_vector;
        concated_mf_int_vector.width(8);
        concated_mf_int_vector.resize(concated_mf.size());
        for (size_t i = 0; i < concated_mf.size(); i++) { concated_mf_int_vector[i] = concated_mf[i]; }
        sdsl::append_zero_symbol(concated_mf_int_vector);
        symbol_table_ = decltype(symbol_table_)(concated_mf_int_vector);
        for (size_t i = 0; i < concated_mf_int_vector.size(); i++) { concated_mf_int_vector[i] = symbol_table_[concated_mf_int_vector[i]]; }
        xbwt->insert(concated_mf_int_vector, symbol_table_.GetEffectiveAlphabetWidth(), symbol_table_.GetAlphabetSize());
    }
}

void tau_lambda_index::load_min_factors(std::string &mf_path) {
    std::ifstream in(mf_path);
    if (!in.is_open()) {
        std::cerr << "Error: Could not open input minimal_factors file " << mf_path << std::endl;
        return;
    }

    std::string delimiters_tmp;
    size_t n;
    in >> inputTextPath >> tau_l >> tau_u >> lambda >> delimiters_tmp >> n;
    // for (auto c : tmp) { delimiters.insert(c); }
    while (n > 0) {
        size_t a, b;
        in >> a >> b;
        min_factors.push_back({a, b});
        n--;
    }

    if (tau_u == 0) { is_original_index = true; }

    if (is_original_index) {
        if (index_type == index_types::old_tau_lambda_type) { throw std::invalid_argument("old tau-lambda-index has no original self-index version"); }
        else if (index_type == index_types::compact_suffix_trie) { throw std::invalid_argument("compact_suffix_trie has no original self-index version"); }
    }

    in.close();
}

// Constructor for generating query patterns
tau_lambda_index::tau_lambda_index(std::string &mf_path) {
    load_min_factors(mf_path);
    load_Text(text, inputTextPath);
}

// Constructor for r-index, LMS and old-tau-lambda-index(DCC version)
tau_lambda_index::tau_lambda_index(std::string &mf_path, index_types index_type): index_type(index_type) {
    load_min_factors(mf_path);
    load_Text(text, inputTextPath);

    if (index_type == index_types::old_tau_lambda_type) {
        old_tau_lambda = new Index(text, min_factors, tau_l, tau_u, lambda);
    } else if (index_type == index_types::compact_suffix_trie) {
        build_XBWT(text);
        location_trie = new compact_suffix_trie(symbol_table_, xbwt.get(), min_factors, text, lambda);
    } else {
        build_XBWT(text);
        std::string maskedTextPath = "tmpMaskedText", maskedText;
        gen_masked_text(text, maskedText);
        if (index_type == index_types::LMS_type) { maskedText += '\0'; }
        output_Text(maskedText, maskedTextPath);

        // generate the self-index
        if (index_type == index_types::r_index_type) {
            std::string input;
            load_Text(input, maskedTextPath);

            r_index = new ri::r_index<>(input, true);
        } else if (index_type == index_types::LMS_type) {
            std::string LMSTemp = "/tmp"; // same as the default value, don't change
            lms = lpg_index(maskedTextPath, LMSTemp, 1, 0.5);
        }

        if (std::remove(maskedTextPath.c_str()) != 0) {
            std::cerr << "Not able to delete the masked text tmp file" << std::endl;
            return;
        }
    }
}

// Constructor for LZ77
tau_lambda_index::tau_lambda_index(std::string &mf_path, std::string &index_path, index_types index_type): index_type(index_type) {
    load_min_factors(mf_path);
    load_Text(text, inputTextPath);

    build_XBWT(text);
    std::string maskedTextPath = "tmpMaskedText", maskedText;
    gen_masked_text(text, maskedText);
    output_Text(maskedText, maskedTextPath);

    if (index_type == index_types::lz77_type) {
        unsigned char br=0;
        unsigned char bs=0;
        unsigned char ss=0;
        char *in = new char[maskedTextPath.size() + 1], *out = new char[index_path.size() + 1];
        std::strcpy(in, maskedTextPath.c_str());
        std::strcpy(out, index_path.c_str());
        lz77 = lz77index::static_selfindex_lz77::build(in, out, br, bs, ss);
        delete[] in, out;
    }

    if (std::remove(maskedTextPath.c_str()) != 0) {
        std::cerr << "Not able to delete the masked text tmp file" << std::endl;
        return;
    }
}

// Constructor for hybrid
tau_lambda_index::tau_lambda_index(std::string &mf_path, std::string &index_path, std::string &outputFolder, index_types index_type, uint M): index_type(index_type) {
    load_min_factors(mf_path);
    load_Text(text, inputTextPath);

    build_XBWT(text);
    std::string maskedTextPath = "tmpMaskedText", maskedText;
    gen_masked_text(text, maskedText);
    output_Text(maskedText, maskedTextPath);

    if (index_type == index_types::hybrid) {
        char *file = new char[maskedTextPath.size()], *out = new char[outputFolder.size()];
        strcpy(file, maskedTextPath.c_str());
        strcpy(out, outputFolder.c_str());
        hybrid_index = new hybrid();
        hybrid_index->build(file, M, out);
        delete[] file, out;
    }

    if (std::remove(maskedTextPath.c_str()) != 0) {
        std::cerr << "Not able to delete the masked text tmp file" << std::endl;
        return;
    }
}

void tau_lambda_index::serialize(std::ofstream &out) {
    size_t length = inputTextPath.size();
    out.write(reinterpret_cast<char*>(&length), sizeof(length));
    out.write(inputTextPath.data(), length);

    sdsl::write_member(index_type, out);
    sdsl::write_member(tau_l, out);
    sdsl::write_member(tau_u, out);
    sdsl::write_member(lambda, out);
    if (index_type == index_types::old_tau_lambda_type) {
        old_tau_lambda->serialize(out);
    } else {
        sdsl::write_member(masked_ratio, out);
        if (tau_u > 0) {
            symbol_table_.Serialize(out);
            xbwt->Serialize(out);
        }
        if (index_type == index_types::r_index_type) {
            r_index->serialize(out);
        } else if (index_type == index_types::LMS_type) {
            lms.serialize(out, NULL, "");
        } else if (index_type == index_types::compact_suffix_trie) {
            location_trie->serialize(out);
        } else if (index_type == index_types::hybrid) {
            hybrid_index->serialized();
        }
    }
}

void tau_lambda_index::load(std::ifstream &in, std::string inputIndexPath, bool xbwt_only) {
    size_t length;
    in.read(reinterpret_cast<char*>(&length), sizeof(length));
    inputTextPath.resize(length);
    in.read(&inputTextPath[0], length);
    sdsl::read_member(index_type, in);
    sdsl::read_member(tau_l, in);
    sdsl::read_member(tau_u, in);
    if (tau_u == 0) { is_original_index = true; }
    sdsl::read_member(lambda, in);
    if (xbwt_only && (index_type == index_types::old_tau_lambda_type || is_original_index)) {
        throw std::invalid_argument("no xbwt for this configuration.");
    }
    if (index_type == index_types::old_tau_lambda_type) {
        load_Text(text, inputTextPath);
        old_tau_lambda = new Index(text);
        old_tau_lambda->load(in);
    } else {
        sdsl::read_member(masked_ratio, in);
        if (!is_original_index) {
            symbol_table_.Load(in);
            xbwt->Load(in);
        }
        if (!xbwt_only) {
            if (index_type == index_types::r_index_type) {
                r_index = new ri::r_index<>();
                r_index->load(in);
            } else if (index_type == index_types::lz77_type) {
                std::string lz77Path = inputIndexPath + "_lz77";
                FILE* fd = fopen(lz77Path.c_str(), "r");
                lz77 = lz77index::static_selfindex_lz77::load(fd);
                fclose(fd);
            } else if (index_type == index_types::LMS_type) {
                lms.load(in);
            } else if (index_type == index_types::compact_suffix_trie) {
                location_trie = new compact_suffix_trie();
                location_trie->load(in);
            } else if (index_type == index_types::hybrid) {
                hybrid_index = new hybrid();
                std::string indexOutFolderDir= "./" + inputIndexPath + "_hybrid/";
                hybrid_index->load(indexOutFolderDir);
            }
        }
    }
}

void tau_lambda_index::_locate(std::string &pattern, std::vector<uint64_t> &results, bool xbwt_only) {
    results.clear();
    sdsl::int_vector<> pattern_int;
    pattern_int.width(8);
    pattern_int.resize(pattern.size());
    size_t p_len = pattern.length();
    for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = symbol_table_[static_cast<unsigned char>(pattern[i]) + t_symbol]; }
    if (xbwt_only) {
        xbwt->match_if_exist(pattern_int);
    } else if (index_type == index_types::old_tau_lambda_type && p_len <= lambda) {
        std::vector<size_t> tmp = old_tau_lambda->location_tree_search(pattern);
        for (auto r : tmp) { results.push_back(r); }
    } else if (p_len <= lambda && index_type == index_types::compact_suffix_trie) {
        size_t offset, length, rank;
        xbwt->match_pos_in_pattern(pattern_int, offset, length, rank);
        if (rank > 0) {
            location_trie->locate(pattern_int, offset, length, rank - 1, tau_l, results); // -1 for shifting to 0-index
        }
    } else if (p_len <= lambda && xbwt->match_if_exist(pattern_int)) {
        if (index_type == index_types::r_index_type) {
            results = r_index->locate_all_tau(pattern, tau_l, tau_u);
        } else if (index_type == index_types::lz77_type) {
            unsigned char *p = new unsigned char[pattern.size() + 1];
            std::memcpy(p, pattern.c_str(), pattern.size() + 1);
            unsigned int nooc;
            std::vector<unsigned int> *tmp = lz77->locate(p, pattern.size(), &nooc);
            if (tmp->size() >= tau_l) {
                results.resize(tmp->size());
                for (size_t i = 0; i < tmp->size(); i++) { results[i] = static_cast<uint64_t>((*tmp)[i]); }
            }
            delete[] p;
        } else if (index_type == index_types::LMS_type) {
            std::set<lpg_index::size_type> tmp;
            lms.locate(pattern, tmp);
            if (tmp.size() >= tau_l) {
                for (auto r : tmp) { results.push_back(r); }
            }
        } else if (index_type == index_types::hybrid) {
            ulong nOcc, *occ;
            uchar* p = new uchar[pattern.size() + 1];
            for (size_t i = 0, e = pattern.size(); i < e; i++) { p[i] = (uchar)(pattern[i]); }
            p[pattern.size()] = '\0';
            hybrid_index->locate(p, lambda, &nOcc, &occ, tau_u);
            if (tau_l <= nOcc) {
                results.resize(nOcc);
                for (size_t i = 0; i < nOcc; i++) { results[i] = occ[i]; }
            }
            delete[] p, occ;
        }
    }
}

void tau_lambda_index::_locate_original_index(std::string &pattern, std::vector<uint64_t> &results, bool xbwt_only) {
    results.clear();
    
    if (index_type == index_types::r_index_type) {
        results = r_index->locate_all_tau(pattern, tau_l, tau_u);
    } else if (index_type == index_types::lz77_type) {
        unsigned char *p = new unsigned char[pattern.size() + 1]; 
        std::memcpy(p, pattern.c_str(), pattern.size() + 1);
        unsigned int nooc;
        std::vector<unsigned int> *tmp = lz77->locate(p, pattern.size(), &nooc, tau_u + 1);
        if (tau_l <= tmp->size() && tmp->size() <= tau_u) {
            results.resize(tmp->size());
            for (size_t i = 0; i < tmp->size(); i++) { results[i] = static_cast<uint64_t>((*tmp)[i]); }
        }
    } else if (index_type == index_types::LMS_type) {
        std::set<lpg_index::size_type> tmp;
        lms.locate(pattern, tmp, tau_u);
        if (tau_l <= tmp.size() && tmp.size() <= tau_u) {
            for (auto r : tmp) { results.push_back(r); }
        }
    } else if (index_type == index_types::hybrid) {
        ulong nOcc, *occ;
        uchar* p = new uchar[pattern.size() + 1];
        for (size_t i = 0, e = pattern.size(); i < e; i++) { p[i] = (uchar)(pattern[i]); }
        p[pattern.size()] = '\0';
        hybrid_index->locate(p, lambda, &nOcc, &occ);
        if (tau_l <= nOcc && nOcc <= tau_u) {
            results.resize(nOcc);
            for (size_t i = 0; i < nOcc; i++) { results[i] = occ[i]; }
        }
        delete[] p;
    }
}

void tau_lambda_index::locate(std::ifstream &in, std::ofstream &out, bool xbwt_only) {
    std::chrono::steady_clock::time_point t1, t2;

    std::string header;
	std::getline(in, header);

    size_t n = get_pattern_info("number=", header);
    size_t m = get_pattern_info("length=", header);
    if (is_original_index) {
        tau_l = get_pattern_info("tau_l=", header);
        tau_u = get_pattern_info("tau_u=", header);
        lambda = get_pattern_info("lambda=", header);
    } else if (tau_l != get_pattern_info("tau_l=", header) ||
               tau_u != get_pattern_info("tau_u=", header) ||
               lambda != get_pattern_info("lambda=", header)) {
        throw std::invalid_argument("tau_l, tau_u, lambda setting mismatch between the index and query patterns");
    }
    size_t total_time = 0, total_cnt = 0;

    for (size_t i = 0; i < n; i++) {
        std::string pattern;
		for(size_t j = 0; j < m; j++){
            char c;
			in.get(c);
			pattern += c;
		}

        std::vector<uint64_t> results;
        if (is_original_index) {
            t1 = std::chrono::steady_clock::now();
            _locate_original_index(pattern, results, xbwt_only);
            t2 = std::chrono::steady_clock::now();
        } else {
            t1 = std::chrono::steady_clock::now();
            _locate(pattern, results, xbwt_only);
            t2 = std::chrono::steady_clock::now();
        }
        
        std::sort(results.begin(), results.end());
        out << "pattern [" << (i+1) << "] -- " << results.size() <<  " occurrences\n";
        if (results.size()) {
            for (auto r : results) { out << r << "\t"; }
            out << "\n";
        }
        out << "time consuming (us): " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "\n\n";
        total_time += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        total_cnt += results.size();
    }

    out << "time(us): " << total_time << "\n";
    out << "occ: " << total_cnt << "\n";
}