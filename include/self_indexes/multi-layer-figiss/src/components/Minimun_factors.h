#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <climits>
#include <string_view>
#include <sdsl/suffix_arrays.hpp>
#include "Alpha_manager.h"

bool debug = false;

const std::unordered_set<char> delemited_chars = {'$', '#'}; // TODO more general way to handle this // TODO move to Alpha_manager
// constexpr unsigned int max_pattern_length = 1000; // manage this as template parameter
// constexpr unsigned int max_pattern_length = 6; // manage this as template parameter

class Minimun_factors {
public:

    // construct option: 0 for double fm, 1 for hash, 2 for lossy hash
    Minimun_factors(const std::string& text, int tau, int lambda = -1, bool activate_lambda = false, int construction_option = 0);
    ~Minimun_factors();

    std::vector<std::pair<size_t, size_t>> get_min_factors();

private:
    std::vector<std::pair<size_t, size_t>> minimun_factors;

    void construct_with_double_fm(const std::string& text, int tau, int lambda, bool activate_lambda);
    void hash_construction(const std::string& text, int tau, int lambda, bool activate_lambda);
    void lossy_hash_construction(const std::string& text, int tau, int lambda, bool activate_lambda);
    void lambda_filter(int lambda, bool activate_lambda);
};

Minimun_factors::Minimun_factors(const std::string& text, int tau, int lambda, bool activate_lambda, int construction_option) {
    if (activate_lambda && lambda < 0) {
        std::cout << "Error: lambda is not set" << std::endl;
        throw std::runtime_error("lambda is not set"); // TODO more specific exception?
    }

    // TODO more general way to handle this
    if (activate_lambda == false) {
        lambda = text.size();
    }

    if (construction_option == 0) {
        construct_with_double_fm(text, tau, lambda, activate_lambda);
    } else if (construction_option == 1) {
        hash_construction(text, tau, lambda, activate_lambda);
    } else if (construction_option == 2) {
        lossy_hash_construction(text, tau, lambda, activate_lambda);
    } else {
        std::cout << "Error: construction option is not supported" << std::endl;
        throw std::runtime_error("construction option is not supported"); // TODO more specific exception?
    }
    lambda_filter(lambda, activate_lambda);
}

Minimun_factors::~Minimun_factors() {}

void Minimun_factors::construct_with_double_fm(const std::string& text, int tau, int lambda, bool activate_lambda) {
    // TODO the meaning of the parameters?
    using fm_index_type = sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, 512, 1024> ;
    using fm_size_type = typename fm_index_type::size_type;

    fm_index_type forward_fm, backward_fm;
    sdsl::construct_im(forward_fm, std::string(text.rbegin(), text.rend()), 1);
    sdsl::construct_im(backward_fm, text, 1);

    // boundary of the interval
    fm_size_type lb = 0, rb = forward_fm.size() - 1;
    for (size_t l = 0, r = 0; r < text.size(); r++) {
        if (delemited_chars.find(text[r]) != delemited_chars.end()) {
            l = r + 1;
            lb = 0;
            rb = forward_fm.size() - 1;
            continue;
        }
        
        auto cnt = sdsl::backward_search(forward_fm, lb, rb, text[r], lb, rb);

        if (cnt <= tau) {
            fm_size_type sub_lb = 0, sub_rb = backward_fm.size() - 1;
            for (size_t i = r; i >= l; --i){
                auto subCnt = sdsl::backward_search(backward_fm, sub_lb, sub_rb, text[i], sub_lb, sub_rb);
                if (subCnt <= tau) {
                    minimun_factors.emplace_back(i, r);
                    l = i + 1;
                    r = i;
                    break;
                }
            }
            lb = 0;
            rb = forward_fm.size() - 1;
        }
    }
}

// rolling hash on string of length lambda
void Minimun_factors::hash_construction(const std::string& text, int tau, int lambda, bool activate_lambda) {
    using hash_key = std::bitset<max_pattern_length * 3>;
    std::string_view text_view(text);
    std::unordered_map<hash_key, std::pair<size_t, size_t>> table; // pair<0>: frequency, pair<1>: length
    
    // sliding window with size lambda
    for (size_t window_begin = 0; window_begin < text.size(); ++window_begin) {
        size_t window_len = (window_begin + lambda < text.size()) ? lambda : text.size() - window_begin;
        std::string_view window(text_view.substr(window_begin, window_len));
        hash_key key = 0;
        
        // insert all prefix in this window
        for (size_t i = 0; i < window.size(); ++i) {
            // if meet a delemited_char, break, TODO more general way to handle this
            if (delemited_chars.find(window[i]) != delemited_chars.end()) {
                break;
            }

            key <<= 3;
            key |= Alpha_manager::char_to_int(window[i]); 

            if (table.find(key) == table.end()) {
                table[key] = {1, i + 1};
            } else {
                table[key].first++;
            }
        }
    }

    std::vector<hash_key> masks;
    hash_key mask("111");
    for (size_t i = 0; i < lambda; ++i) {
        masks.push_back(~(mask << (i * 3)));
    }
    
    // remove factors that are not minimun
    std::unordered_set<hash_key> minimun_factors;
    for (const auto& [key, val] : table) {
        auto [freq, len] = val;
        if (freq <= tau) {
            hash_key prefix = key >> 3;
            hash_key suffix = key & masks[len - 1];

            if (debug) { // TODO debug
                // std::cout << key << std::endl;
                // std::cout << prefix << std::endl;
                // std::cout << suffix << std::endl;
                std::cout << Alpha_manager::hash_key_to_string(key) << " " << freq << " " << len << std::endl;
                std::cout << Alpha_manager::hash_key_to_string(prefix) << " " << table[prefix].first << " " << table[prefix].second << std::endl;
                std::cout << Alpha_manager::hash_key_to_string(suffix) << " " << table[suffix].first << " " << table[suffix].second << std::endl;
                std::cout << std::endl;
            }

            if (len == 1 || (table[prefix].first > tau && table[suffix].first > tau)) {
                minimun_factors.insert(key);
            }
        }
    }
    
    if (debug) { // TODO debug
        for (auto key : minimun_factors) {
            // back to string
            std::cout << Alpha_manager::hash_key_to_string(key) << std::endl;
        }
    }

    // slide the window again to parse all minimun factors in text
    std::vector<std::pair<size_t, size_t>> minimun_factors_index;
    for (size_t window_begin = 0; window_begin < text.size(); ++window_begin) {
        size_t window_len = (window_begin + lambda < text.size()) ? lambda : text.size() - window_begin;
        std::string_view window(text_view.substr(window_begin, window_len));
        hash_key key;
        
        // check all prefix in this window
        for (size_t i = 0; i < window.size(); ++i) {
            if (delemited_chars.find(window[i]) != delemited_chars.end()) {
                break;
            }

            key <<= 3;
            key |= Alpha_manager::char_to_int(window[i]);

            if (minimun_factors.find(key) != minimun_factors.end()) {
                minimun_factors_index.emplace_back(window_begin, window_begin + i);
            }
        }
    }

    this->minimun_factors = minimun_factors_index;
}

void Minimun_factors::lossy_hash_construction(const std::string& text, int tau, int lambda, bool activate_lambda) {
    using hash_key_type = unsigned long long;
    using table_type = std::unordered_map<hash_key_type, std::vector<std::pair<size_t, size_t>>>;

    table_type prev, curr;
    std::unordered_set<hash_key_type> prev_few_than_tau;
    std::vector<std::pair<size_t, size_t>> minimun_factors_index;

    hash_key_type base = 7;
    // unsigned long mod = 4294967291; // TODO better mod?
    hash_key_type mod = LLONG_MAX / 10;

    // cases for len = 1
    for (size_t i = 0; i < text.size(); ++i) {
        if (delemited_chars.find(text[i]) != delemited_chars.end()) {
            continue;
        }
        hash_key_type key = Alpha_manager::char_to_int(text[i]);
        if (curr.find(key) == curr.end()) {
            curr[key] = {{i, i}};
        } else {
            curr[key].emplace_back(i, i);
        }
    }
    for (auto& [key, vec] : curr) {
        if (vec.size() <= tau) {
            minimun_factors_index.insert(minimun_factors_index.end(), vec.begin(), vec.end());
            prev_few_than_tau.insert(key);
            continue;
        }
        prev[key] = vec; // TODO move?
    }

    // cases for len > 1
    for (size_t len = 2; len <= lambda; ++len) {
        curr.clear();
        for (const auto& [key, vec] : prev) {
            for (auto [l, r] : vec) {
                if (r + 1 >= text.size()) continue;
                if (delemited_chars.find(text[r + 1]) != delemited_chars.end()) {
                    continue;
                }
                hash_key_type new_key = (key * base + Alpha_manager::char_to_int(text[r + 1])) % mod;
                if (curr.find(new_key) == curr.end()) {
                    curr[new_key] = {{l, r + 1}};
                } else {
                    curr[new_key].emplace_back(l, r + 1);
                }
            }
        }

        // calculate base_pow by direct calculation to prevent floating point error
        hash_key_type base_pow = 1;
        for (int i = 0; i < len - 1; ++i) {
            base_pow *= base;
            base_pow %= mod;
        }

        std::unordered_set<hash_key_type> curr_few_than_tau;
        prev.clear();
        for (const auto& [key, vec] : curr) {
            if (vec.size() > tau) {
                prev[key] = vec; // TODO move?
                continue;
            }

            // few than tau
            size_t left = vec[0].first;
            hash_key_type most_significant_digit = (Alpha_manager::char_to_int(text[left]) * base_pow) % mod;
            hash_key_type suffix_key = (key + mod - most_significant_digit) % mod;
            if (prev_few_than_tau.find(suffix_key) == prev_few_than_tau.end()) {
                minimun_factors_index.insert(minimun_factors_index.end(), vec.begin(), vec.end());
            }
            curr_few_than_tau.insert(key);
        }
        prev_few_than_tau = curr_few_than_tau;
    }

    this->minimun_factors = minimun_factors_index;
}


// remove factors that are longer than lambda
// TODO inline to the constructor?
void Minimun_factors::lambda_filter(int lambda, bool activate_lambda) {
    if (activate_lambda) {
        auto it = remove_if(minimun_factors.begin(), minimun_factors.end(), 
            [lambda](const std::pair<size_t, size_t>& p) {return p.second - p.first + 1 > lambda;});
        minimun_factors.erase(it, minimun_factors.end());
    }
}

std::vector<std::pair<size_t, size_t>> Minimun_factors::get_min_factors() {
    return minimun_factors;
}