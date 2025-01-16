// manage the conversion between char, int, string, and hash_key(bitset)
#pragma once
#include <bitset>
#include <unordered_map>

constexpr int max_pattern_length = 10; // TODO to template
using hash_key = std::bitset<max_pattern_length * 3>; // TODO to template


// std::unordered_map<char, int> char_to_int_map = { // TODO for debug
//     {'$', 1},
//     {'S', 2}, {'R', 3}, {'D', 4}, {'C', 5}, {'B', 6}, {'A', 7},
// };

// std::unordered_map<int, char> int_to_char_map = { // TODO for debug
//     {1, '$'},
//     {2, 'S'}, {3, 'R'}, {4, 'D'}, {5, 'C'}, {6, 'B'}, {7, 'A'},
// };

std::unordered_map<char, int> char_to_int_map = { // TODO unsigned int?
    {'$', 0}, {'#', 1}, 
    {'A', 2}, {'C', 3}, {'G', 4}, {'T', 5}, {'N', 6}, 
    {'a', 2}, {'c', 3}, {'g', 4}, {'t', 5}, {'n', 6}
};

std::unordered_map<int, char> int_to_char_map = { // TODO unsigned int?
    {0, '$'}, {1, '#'}, 
    {2, 'A'}, {3, 'C'}, {4, 'G'}, {5, 'T'}, {6, 'N'}
};

class Alpha_manager {
private:
    // static std::unordered_map<char, int> char_to_int_map =  {
    //     {'$', 0}, {'#', 1}, 
    //     {'A', 2}, {'C', 3}, {'G', 4}, {'T', 5}, {'N', 6}, 
    //     {'a', 2}, {'c', 3}, {'g', 4}, {'t', 5}, {'n', 6}
    // };
    // static std::unordered_map<int, char> int_to_char_map = {
    //     {0, '$'}, {1, '#'}, 
    //     {2, 'A'}, {3, 'C'}, {4, 'G'}, {5, 'T'}, {6, 'N'}
    // };

public:
    static std::string hash_key_to_string(const hash_key& key) {
        hash_key temp = key;
        std::string res;
        while (temp.any()) {
            int curr_char = (temp & hash_key("111")).to_ulong();
            temp = temp >> 3;
            res += int_to_char(curr_char);
        }
        return res;
    }

    // TODO debug
    // static hash_key string_to_hash_key(const std::string& str) {
    //     hash_key result;
    //     for (auto& c : str) {
    //         result.set(char_to_int.at(c));
    //     }
    //     return result;
    // }

    static int char_to_int(char c) {
        return char_to_int_map.at(c);
    }

    static char int_to_char(int i) {
        return int_to_char_map.at(i); // exception when i is not in the map
    }

    static void set_char_to_int_map(const std::unordered_map<char, int>& map) {
        char_to_int_map = map;
    }

    static void set_int_to_char_map(const std::unordered_map<int, char>& map) {
        int_to_char_map = map;
    }

    static void reset_char_to_int_map() {
        char_to_int_map = {
            {'$', 0}, {'#', 1}, 
            {'A', 2}, {'C', 3}, {'G', 4}, {'T', 5}, {'N', 6}, 
            {'a', 2}, {'c', 3}, {'g', 4}, {'t', 5}, {'n', 6}
        };
    }

    static void reset_int_to_char_map() {
        int_to_char_map = {
            {0, '$'}, {1, '#'}, 
            {2, 'A'}, {3, 'C'}, {4, 'G'}, {5, 'T'}, {6, 'N'}
        };
    }
};