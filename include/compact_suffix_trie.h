#pragma once

#include <string>
#include <vector>
#include <sdsl/int_vector.hpp>
#include "xbwt/xbwt.h"
#include "xbwt/xbwt_location.h"
#include "symbol_table/symbol_table.h"

class compact_suffix_trie {
private:
    // XBWT_location *xbwt_forw {new XBWT_location()}, *xbwt_reverse {new XBWT_location()};
    uint64_t t_symbol = static_cast<uint64_t>(1); // terminate symbols
    sdsl::int_vector<> AC_auto_leaves_offsets;
    std::vector<XBWT_location*> xbwts_forw, xbwts_reverse;
    // sdsl::int_vector<> xbwt_forw_locations, xbwt_forw_leaf_offsets, xbwt_reverse_locations, xbwt_reverse_leaf_offsets;

    void build_forw_reverse_xbwts(SymbolTable &symbol_table_, XBWT* xbwt, std::vector<std::pair<size_t, size_t>> &min_factors, std::string &text, size_t lambda);
    void gen_local_locations(std::vector<unsigned char> &inserted_string, std::vector<size_t> &inserted_locations, std::vector<size_t> &leaves_locations, XBWT_location *xbwt_used);

public:
    compact_suffix_trie(SymbolTable &symbol_table_, XBWT* xbwt, std::vector<std::pair<size_t, size_t>> &min_factors, std::string &text, size_t lambda);
    void serialize(std::ofstream &out) {
        AC_auto_leaves_offsets.serialize(out);
        sdsl::write_member(xbwts_forw.size(), out);
        for (auto xbwt : xbwts_forw) { xbwt->Serialize(out); }
        for (auto xbwt : xbwts_reverse) { xbwt->Serialize(out); }
    }
};

compact_suffix_trie::compact_suffix_trie(SymbolTable &symbol_table_, XBWT* xbwt, std::vector<std::pair<size_t, size_t>> &min_factors, std::string &text, size_t lambda) {
    build_forw_reverse_xbwts(symbol_table_, xbwt, min_factors, text, lambda);
}

void compact_suffix_trie::build_forw_reverse_xbwts(SymbolTable &symbol_table_, XBWT* AC_auto, std::vector<std::pair<size_t, size_t>> &min_factors, std::string &text, size_t lambda) {
    size_t n = text.size(), leaf_cnt = AC_auto->getGrammarNumber();
    AC_auto_leaves_offsets.resize(leaf_cnt);
    xbwts_forw = std::vector<XBWT_location*>(leaf_cnt, new XBWT_location());
    xbwts_reverse = std::vector<XBWT_location*>(leaf_cnt, new XBWT_location());
    
    // concatenated extended substrings and collect the locations (recording b of each minimal factor, [b, e])
    std::vector<std::vector<unsigned char>> strings_forw(leaf_cnt), strings_reverse(leaf_cnt);
    std::vector<std::vector<size_t>> locations_tmp(leaf_cnt);
    
    for (auto [b, e] : min_factors) {
        sdsl::int_vector<> pattern_int;
        pattern_int.width(8);
        pattern_int.resize(e + 1 - b);
        for (size_t i = b; i <= e; i++) { pattern_int[i] = symbol_table_[static_cast<unsigned char>(text[i]) + t_symbol]; }
        size_t rank = AC_auto->match(pattern_int.begin(), pattern_int.end()) - 1; // shift to 0-index

        size_t forw_l = 0;
        if (e + 1 > lambda) { forw_l = e - lambda + 1; }
        for (size_t i = forw_l; i < b; i++) {
            strings_forw[rank].push_back(symbol_table_[static_cast<unsigned char>(text[i]) + t_symbol]);
        }
        strings_forw[rank].push_back(symbol_table_[t_symbol]);
        size_t reverse_r = n - 1;
        if (b + lambda < n + 1) { reverse_r = b + lambda - 1; }
        for (size_t i = reverse_r; i > e; i--) {
            strings_reverse[rank].push_back(symbol_table_[static_cast<unsigned char>(text[i]) + t_symbol]);
        }
        strings_reverse[rank].push_back(symbol_table_[t_symbol]);

        locations_tmp[rank].push_back(b);
    }

    // build xbwts for each leaf in AC_auto
    for (size_t i = 0; i < leaf_cnt; i++) {
        sdsl::int_vector<> forw_int_vector, reverse_int_vector;
        forw_int_vector.width(8); // TODO: change to the smallest bits usage
        forw_int_vector.resize(strings_forw[i].size());
        for (size_t j = 0; j < strings_forw[i].size(); j++) { forw_int_vector[i] = strings_forw[i][j]; }
        sdsl::append_zero_symbol(forw_int_vector);
        xbwts_forw[i]->insert(forw_int_vector, symbol_table_.GetEffectiveAlphabetWidth(), symbol_table_.GetAlphabetSize());
        reverse_int_vector.width(8); // TODO: change to the smallest bits usage
        reverse_int_vector.resize(strings_reverse[i].size());
        for (size_t j = 0; j < strings_reverse[i].size(); j++) { reverse_int_vector[i] = strings_reverse[i][j]; }
        sdsl::append_zero_symbol(reverse_int_vector);
        xbwts_reverse[i]->insert(reverse_int_vector, symbol_table_.GetEffectiveAlphabetWidth(), symbol_table_.GetAlphabetSize());
    }

    // locations corresponding to each xbwt
    //  [b_forw, e_forw) is the text to be matched (text[e_forw] == t_symbol)
    for (size_t i = 0, cnt = 0, b_forw = 0, e_forw = 0, b_rev = 0, transformed_t_symbol = symbol_table_[t_symbol], e_rev = 0; i < leaf_cnt; i++) {
        AC_auto_leaves_offsets[i] = cnt;
        cnt += locations_tmp[i].size();

        // std::vector<size_t> leaves_locations;
        size_t xbwt_leaves_cnt = xbwts_forw[i]->getGrammarNumber();
        std::vector<std::vector<size_t>> xbwt_forw_locations_tmp (xbwt_leaves_cnt), xbwt_reverse_locations_tmp (xbwt_leaves_cnt);;
        std::vector<size_t> xbwt_forw_leaf_offsets_tmp (xbwt_leaves_cnt), xbwt_reverse_leaf_offsets_tmp (xbwt_leaves_cnt);
        size_t idx_forw = 0, idx_rev = 0;
        
        // forward
        while (e_forw < n) {
            while (strings_forw[i][e_forw] != transformed_t_symbol) { e_forw++; }
            sdsl::int_vector<> pattern_forw;
            pattern_forw.width(8);
            pattern_forw.resize(e_forw - b_forw + 1);
            for (size_t j = b_forw, k = pattern_forw.size(); j < e_forw; j++, k--) {
                pattern_forw[k] = strings_forw[i][j];
            }
            size_t rank = xbwts_forw[i]->match(pattern_forw.begin(), pattern_forw.end()) - 1; // shift to 0-index
            xbwt_forw_locations_tmp[rank].push_back(locations_tmp[i][idx_forw]);
            idx_forw++;
            b_forw = ++e_forw;
        }
        
        for (size_t j = 0, cnt = 0; j < xbwt_leaves_cnt; j++) {
            xbwt_forw_leaf_offsets_tmp[j] = cnt;
            cnt += xbwt_forw_locations_tmp[j].size();
        }
        xbwts_forw[i]->insert_locations(xbwt_forw_locations_tmp, locations_tmp.size());

        // reverse
        while (e_rev < n) {
            while (strings_reverse[i][e_rev] != transformed_t_symbol) { e_rev++; }
            sdsl::int_vector<> pattern_reverse;
            pattern_reverse.width(8);
            pattern_reverse.resize(e_rev - b_rev + 1);
            for (size_t j = b_rev, k = pattern_reverse.size(); j < e_rev; j++, k--) {
                pattern_reverse[k] = strings_reverse[i][j];
            }
            size_t rank = xbwts_reverse[i]->match(pattern_reverse.begin(), pattern_reverse.end()) - 1; // shift to 0-index
            xbwt_reverse_locations_tmp[rank].push_back(locations_tmp[i][idx_rev]);
            idx_rev++;
            b_rev = ++e_rev;
        }
        
        for (size_t j = 0, cnt = 0; j < xbwt_leaves_cnt; j++) {
            xbwt_reverse_leaf_offsets_tmp[j] = cnt;
            cnt += xbwt_reverse_locations_tmp[j].size();
        }
    }
}

void compact_suffix_trie::gen_local_locations(std::vector<unsigned char> &inserted_string, std::vector<size_t> &inserted_locations, std::vector<size_t> &leaves_locations, XBWT_location *xbwt_used) {
    leaves_locations.clear();
    
    if (inserted_string.size() == 1) { // only 1 t_symbol
        leaves_locations.push_back(inserted_locations[0]);
    } else {
        // inserted in the reverse order, so we need to query back
        size_t b = inserted_string.size() - 2, e = b;
        for (size_t i = inserted_locations.size(); i > 0; i--) {
            while (e > 0 && inserted_string[e] != t_symbol) { e--; }
            sdsl::int_vector<> pattern_int;
            pattern_int.width(8);
            pattern_int.resize(b - e + 1);
        }

    }
}