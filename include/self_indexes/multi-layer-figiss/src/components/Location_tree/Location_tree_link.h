#pragma once
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "Location_tree_node.h"

// decalare Location_tree_node externally
// TODO better way to do this?
struct Location_tree_node;

struct Location_tree_link {
    uint32_t start;
    uint32_t end;
    std::unique_ptr<Location_tree_node> node;
    Location_tree_link(size_t start, size_t end, size_t factor_index) 
        : start(start), end(end), node(std::make_unique<Location_tree_node>(factor_index)) {};
    Location_tree_link(size_t start, size_t end, std::unique_ptr<Location_tree_node>&& node) 
        : start(start), end(end), node(std::move(node)) {};
};