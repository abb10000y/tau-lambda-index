#pragma once
#include <memory>
#include <map>
#include "Blind_tree_link.h"

struct Blind_tree_link;

struct Blind_tree_node {
    std::map<char, std::unique_ptr<Blind_tree_link>> links; 
    bool is_leaf = false;
    size_t factor_index;
    // bool is_leaf; // TODO remove this or not?
    Blind_tree_node(size_t factor_index, bool is_leaf) : factor_index(factor_index), is_leaf(is_leaf) {};
    Blind_tree_node(size_t factor_index) : factor_index(factor_index), is_leaf(true) {};
    Blind_tree_node() : is_leaf(true) {};
};