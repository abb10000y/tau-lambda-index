// the top layer of the index
#pragma once

#include "Ac_automata.h"
#include "Minimun_factors.h"
#include "Blind_tree/Blind_tree.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include <string_view>
#include <exception>
#include <filesystem>
#include <unordered_set>
#include <sdsl/io.hpp>

class Index {
public:
    // load the serialized AC + LT
    Index(std::string &text): text(text), forward_location_tree(text, false), reverse_location_tree(text, true){}
    Index(std::string &text, const std::vector<std::pair<size_t, size_t>> &mf, size_t tau_l, size_t tau_u, size_t lambda);
    ~Index();

    std::vector<size_t> location_tree_search(const std::string& pattern); // search the pattern in the text using location tree

    void serialize(std::ofstream &out) {
        sdsl::write_member(tau_l, out);
        sdsl::write_member(tau_u, out);
        sdsl::write_member(lambda, out);
        std::cout << "Serializing AC_automaton" << std::endl;
        ac_automata.serialize(out);
        std::cout << "Serializing forward_location_tree" << std::endl;
        forward_location_tree.serialize(out);
        std::cout << "Serializing reverse_location_tree" << std::endl;
        reverse_location_tree.serialize(out);
    }

    void load(std::ifstream &in) {
        sdsl::read_member(tau_l, in);
        sdsl::read_member(tau_u, in);
        sdsl::read_member(lambda, in);
        ac_automata.load(in);
        forward_location_tree.load(in);
        reverse_location_tree.load(in);
        forward_location_tree.set_text(text);
        reverse_location_tree.set_text(text);
    }

private:
    std::string &text; // the text
    size_t tau_u, tau_l, lambda;
    Ac_automata ac_automata; // the ac_automata
    Blind_tree forward_location_tree; // the forward location tree
    Blind_tree reverse_location_tree; // the reverse location tree
};

Index::Index(std::string &text, const std::vector<std::pair<size_t, size_t>> &mf, size_t tau_l, size_t tau_u, size_t lambda):
text(text), tau_l(tau_l), tau_u(tau_u), lambda(lambda), ac_automata(text, mf), forward_location_tree(text, false), reverse_location_tree(text, true)
{
    std::cout << "Inserting min_factors to forward_location_tree" << std::endl;
    forward_location_tree.insert_factors(mf);
    std::cout << "Inserting min_factors to reverse_location_tree" << std::endl;
    reverse_location_tree.insert_factors(mf);
}

Index::~Index() {}

// Blind tree version
std::vector<size_t> Index::location_tree_search(const std::string& pattern) {   
    auto [start_pattern, length] = ac_automata.match_pos_in_pattern(pattern);
    if (start_pattern == -1) {
        std::cout << start_pattern << ", " << length << std::endl;
        return std::vector<size_t>();
    }
    
    std::string_view pattern_view(pattern); // convert pattern to string_view for better performance
    std::vector<size_t> results;
    
    auto factor_end_in_pattern = start_pattern + length - 1; // end of factor in pattern
    std::string pattern_suffix = pattern.substr(start_pattern, pattern.size() - start_pattern); // begin of factor to the end of pattern
    std::string pattern_prefix = pattern.substr(0, factor_end_in_pattern + 1); // begin of pattern to end of factor TODO +1 or not?
    
    std::vector<size_t> forward_factors = forward_location_tree.match(pattern_suffix, pattern_suffix.size(), start_pattern);
    std::vector<size_t> reverse_factors = reverse_location_tree.match(pattern_prefix, pattern_prefix.size(), start_pattern);
    
    std::cout << start_pattern << ", " << length << ", " << forward_factors.size() << ", " << reverse_factors.size() << std::endl;

    std::sort(forward_factors.begin(), forward_factors.end());
    std::sort(reverse_factors.begin(), reverse_factors.end());
    results.resize(std::min(forward_factors.size(), reverse_factors.size()));
    auto it = std::set_intersection(forward_factors.begin(), forward_factors.end(),
    reverse_factors.begin(), reverse_factors.end(),
    results.begin());
    if (it - results.begin() <tau_l) {
        results.clear();
    } else {
        results.resize(it - results.begin());
    }
    
    return results;
}