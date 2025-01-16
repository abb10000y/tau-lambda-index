#pragma once
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "Location_tree_link.h"

// decalare Location_tree_link externally
// TODO better way to do this?
struct Location_tree_link;

struct Location_tree_node {
    std::map<char, std::unique_ptr<Location_tree_link>> links;
    uint32_t factor_index;
    bool is_leaf;
    Location_tree_node() : is_leaf(false) {};
    Location_tree_node(size_t factor_index) : factor_index(factor_index), is_leaf(true) {};
};