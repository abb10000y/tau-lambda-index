#pragma once
#include <map>
#include <memory>
#include "Blind_tree_node.h"

struct Blind_tree_node;

struct Blind_tree_link {
    size_t len;
    std::unique_ptr<Blind_tree_node> node;
    // Blind_tree_link(size_t len, size_t factor_index) : len(len), node(std::make_unique<Blind_tree_node>(factor_index)) {};
    // Blind_tree_link(size_t len) : len(len), node(std::make_unique<Blind_tree_node>()) {};
    Blind_tree_link(size_t len, Blind_tree_node* node) : len(len), node(node) {};
};