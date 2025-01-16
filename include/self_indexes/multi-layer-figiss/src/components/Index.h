// the top layer of the index
#pragma once

#include "Ac_automata.h"
#include "Minimun_factors.h"
#include "Location_tree/Location_tree.h"
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
#include <chrono>
#include "../components/Load_mf_from_disk.h"
#include <sys/resource.h>
#include <sys/time.h>

class Index {
public:
    // generate min_factor from the scratch
    Index(const std::string& text, int tau, int lambda = -1, bool activate_lambda = false, bool use_location_tree = false);
    // load the serialized AC + LT
    Index(std::string& directory_name, const std::string& text, int tau, int lambda = -1, bool activate_lambda = false, bool use_location_tree = false);
    // load the min_factor file
    Index(const std::string& text, int tau, int lambda, bool activate_lambda, bool use_location_tre, std::string mfPath);
    ~Index();

    std::vector<size_t> search(const std::string& pattern, bool verbose = true); // search the pattern in the text
    std::vector<size_t> direct_search(const std::string& pattern); // search the pattern in the text directly
    std::vector<size_t> location_tree_search(const std::string& pattern); // search the pattern in the text using location tree

    void serialize(const std::string& directory_name) {
        namespace fs = std::filesystem;

        // if no such directory, create one

        if (!fs::exists(directory_name)) {
            fs::create_directory(directory_name);
        }

        std::string ac_filename = directory_name + "/ac_automata.bin";
        std::string forward_location_tree_filename = directory_name + "/forward_location_tree.bin";
        std::string reverse_location_tree_filename = directory_name + "/reverse_location_tree.bin";
        ac_automata.serialize(ac_filename);
        forward_location_tree.serialize(forward_location_tree_filename);
        reverse_location_tree.serialize(reverse_location_tree_filename);
    }

    // prop for debug TODO remove
    size_t get_ac_size() {
        return ac_automata.get_ac_size();
    }

private:
    const std::string& text; // the text
    int tau; // the frequency of the index
    int lambda; // the max length of the pattern
    bool activate_lambda; // whether activate lambda
    bool use_location_tree; // whether use location tree
    Ac_automata ac_automata; // the ac_automata
    // Location_tree forward_location_tree; // the forward location tree
    // Location_tree reverse_location_tree; // the reverse location tree
    Blind_tree forward_location_tree; // the forward location tree
    Blind_tree reverse_location_tree; // the reverse location tree
};

// TODO ugly, need to be refactored
Index::Index(const std::string& text, int tau, int lambda, bool activate_lambda, bool use_location_tree)
    : text(text), 
      tau(tau), 
      lambda(lambda), 
      activate_lambda(activate_lambda),
      use_location_tree(use_location_tree),
      forward_location_tree(text, false), // false/true for forward/reverse
      reverse_location_tree(text, true) 
{
    if (activate_lambda && lambda < 0) {
        std::cout << "Error: lambda is not set" << std::endl;
        throw std::runtime_error("lambda is not set"); // TODO more specific exception?
    
    }

    Minimun_factors min_factors(text, tau, lambda, activate_lambda);
    ac_automata = Ac_automata(text, min_factors.get_min_factors());
    forward_location_tree.insert_factors(min_factors.get_min_factors());
    reverse_location_tree.insert_factors(min_factors.get_min_factors());
}

Index::Index(std::string& directory_name, const std::string& text, int tau, int lambda, bool activate_lambda, bool use_location_tree) 
    : text(text), 
      tau(tau), 
      lambda(lambda), 
      activate_lambda(activate_lambda), 
      use_location_tree(use_location_tree),
      forward_location_tree(text, false), // false/true for forward/reverse
      reverse_location_tree(text, true) 
{
    if (activate_lambda && lambda < 0) {
        std::cout << "Error: lambda is not set" << std::endl;
        throw std::runtime_error("lambda is not set"); // TODO more specific exception?
    
    }

    std::string ac_filename = directory_name + "/ac_automata.bin";
    std::string forward_location_tree_filename = directory_name + "/forward_location_tree.bin";
    std::string reverse_location_tree_filename = directory_name + "/reverse_location_tree.bin";
    ac_automata.load(ac_filename);
    forward_location_tree.load(forward_location_tree_filename);
    reverse_location_tree.load(reverse_location_tree_filename);
}

// load the min_factor file
Index::Index(const std::string& text, int tau, int lambda, bool activate_lambda, bool use_location_tree, std::string mfPath)
    : text(text), 
      tau(tau), 
      lambda(lambda), 
      activate_lambda(activate_lambda), 
      use_location_tree(use_location_tree),
      forward_location_tree(text, false), // false/true for forward/reverse
      reverse_location_tree(text, true) 
{
       
    if (activate_lambda && lambda < 0) {
        std::cout << "Error: lambda is not set" << std::endl;
        throw std::runtime_error("lambda is not set"); // TODO more specific exception?
    
    }   

    Load_mf_from_disk mfSet(mfPath);
    ac_automata = Ac_automata(text, mfSet.get_min_factors());
    forward_location_tree.insert_factors(mfSet.get_min_factors());
    reverse_location_tree.insert_factors(mfSet.get_min_factors());
}

Index::~Index() {}

std::vector<size_t> Index::search(const std::string& pattern, bool verbose) {
    if (activate_lambda && pattern.size() > lambda) {
        if (verbose) {
            std::cout << "Warning: the length of the pattern is larger than lambda" << std::endl;
        }
        return std::vector<size_t>();
    }

    if (use_location_tree) {
        return location_tree_search(pattern);
    } 
    else {
        return direct_search(pattern);
    }
}

std::vector<size_t> Index::direct_search(const std::string& pattern) {
    std::vector<size_t> result;
    std::vector<std::tuple<size_t, size_t, size_t>> ac_result = ac_automata.match(pattern);
    std::string_view pattern_view(pattern); // convert pattern to string_view for better performance

    for (auto [start_text, start_pattern, length] : ac_result) {
        size_t begin_pattern_in_text = start_text - start_pattern;
        std::string_view pattern_in_text = std::string_view(text).substr(begin_pattern_in_text, pattern.size());
        if (pattern_in_text == pattern_view) {
            result.emplace_back(begin_pattern_in_text);
        }
    }
    return result;
}

// // TODO str -> string_view for better performance
// // the whole alg need performance improvement
// std::vector<size_t> Index::location_tree_search(const std::string& pattern) {
//     std::set<size_t> result; // TODO set -> unordered_set for better performance
//     std::vector<std::tuple<size_t, size_t, size_t>> ac_result = ac_automata.match(pattern);
//     std::string_view pattern_view(pattern); // convert pattern to string_view for better performance

//     for (auto [dummy, start_pattern, length] : ac_result) {
//         auto factor_end_in_pattern = start_pattern + length - 1; // end of factor in pattern
//         std::string pattern_suffix = pattern.substr(start_pattern, pattern.size() - start_pattern); // begin of factor to the end of pattern
//         std::string pattern_prefix = pattern.substr(0, factor_end_in_pattern + 1); // begin of pattern to end of factor TODO +1 or not?

//         std::vector<size_t> forward_factors = forward_location_tree.match(pattern_suffix);
//         std::vector<size_t> reverse_factors = reverse_location_tree.match(pattern_prefix);
//         std::set<size_t> forward_result(forward_factors.begin(), forward_factors.end()); // TODO set -> unordered_set for better performance

//         for (auto start_text : reverse_factors) {
//             if (forward_result.find(start_text) != forward_result.end()) {
//                 size_t begin_pattern_in_text = start_text - start_pattern;
//                 result.emplace(begin_pattern_in_text);
//             }
//         }
//     }
//     return std::vector<size_t>(result.begin(), result.end());
// }

// Blind tree version
std::vector<size_t> Index::location_tree_search(const std::string& pattern) {
    // std::vector<std::tuple<size_t, size_t, size_t>> ac_result = ac_automata.match(pattern);

    //std::chrono::steady_clock::time_point t1, t2, t3, t4, t5;
    //std::ofstream outputFile("./output_partitioin/summary_NotExist_with_mf.txt", std::ios::app);
    // outputFile << "AC matching (ns)" << "\t" << "gen LT (ns)" << "\t" << "loop reverse_factors (ns)" << "\t" << "loop candidates_result (ns)" << "\n";

    //t1 = std::chrono::steady_clock::now();
    
    auto [start_pattern, length] = ac_automata.match_pos_in_pattern(pattern);
    if (start_pattern == -1) {
        return std::vector<size_t>();
    }

    //t2 = std::chrono::steady_clock::now();

    std::string_view pattern_view(pattern); // convert pattern to string_view for better performance
    std::unordered_set<size_t> candidate_result;

    auto factor_end_in_pattern = start_pattern + length - 1; // end of factor in pattern
    std::string pattern_suffix = pattern.substr(start_pattern, pattern.size() - start_pattern); // begin of factor to the end of pattern
    std::string pattern_prefix = pattern.substr(0, factor_end_in_pattern + 1); // begin of pattern to end of factor TODO +1 or not?

    std::vector<size_t> forward_factors = forward_location_tree.match(pattern_suffix);
    std::vector<size_t> reverse_factors = reverse_location_tree.match(pattern_prefix);
    std::set<size_t> forward_result(forward_factors.begin(), forward_factors.end()); // TODO set -> unordered_set for better performance

    //t3 = std::chrono::steady_clock::now();

    for (auto start_text : reverse_factors) {
        if (forward_result.find(start_text) != forward_result.end()) {
            size_t begin_pattern_in_text = start_text - start_pattern;
            candidate_result.emplace(begin_pattern_in_text);
        }
    }

    //t4 = std::chrono::steady_clock::now();

    // for each candidate result, check if it is a real result
    std::vector<size_t> result;
    for (auto begin_pattern_in_text: candidate_result) {
        std::string_view pattern_in_text = std::string_view(text).substr(begin_pattern_in_text, pattern.size());
        if (pattern_in_text == pattern_view) {
            result.emplace_back(begin_pattern_in_text);
        }
    }

    //t5 = std::chrono::steady_clock::now();

    //outputFile << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << "\t";
    //outputFile << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() << "\t";
    //outputFile << std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count() << "\t";
    //outputFile << std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count() << "\n";
    //outputFile.close();

    return result;
}