// Compact trie structure for location matching
// TODO combine the forward and reverse version

#pragma once
#include <string>
#include <vector>
#include <stack>
#include <tuple>
#include <utility>
#include <memory>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <functional> // TODO remove
#include "../Location_tree/Location_tree.h"
#include "Blind_tree_node.h"
#include "Blind_tree_link.h"
#include <sdsl/bit_vectors.hpp>

class Blind_tree {
public:
    using Node = Blind_tree_node;
    using Link = Blind_tree_link;

    // construct a location tree from the text
    // reverse: whether to construct a reverse location tree
    Blind_tree(std::string& text, bool reverse = false);
    ~Blind_tree();

    void insert_factors(const std::vector<std::pair<size_t, size_t>>& min_factors);
    // std::vector<size_t> match(const std::string& pattern);
    std::vector<size_t> match(std::string_view pattern, size_t factor_size, size_t offset);

    // debug, print info of the tree for each node
    void print(Node* node = nullptr, bool first_call = true, bool print_out_path = true); 
    size_t get_num_nodes(); // TODO remove

    void serialize(std::ofstream& out);
    void load(std::ifstream& in);
    void set_text(std::string &T) { text = T; }

private:
    std::unique_ptr<Node> root;
    std::string& text;
    bool reverse;
    bool debug = false; // TODO remove

    // recursively get all the locations of the subtree
    // TODO non dfs version
    void get_locations(Node* node, std::vector<size_t>& result, size_t offset); 

    // dfs to pick left most location
    size_t pick_one_location(Node* node);
};

Blind_tree::Blind_tree(std::string& text, bool reverse) : text(text), reverse(reverse) {}

Blind_tree::~Blind_tree() {}

void Blind_tree::insert_factors(const std::vector<std::pair<size_t, size_t>>& min_factors) {
    Location_tree location_tree(text, reverse);
    location_tree.insert_factors(min_factors);

    // dfs function def
    std::function<Node*(Location_tree_node*)> dfs = [&](Location_tree_node* loc_curr) -> Node* {
        if (!loc_curr) {
            return nullptr;
        }
        auto node = std::make_unique<Node>(loc_curr->factor_index, loc_curr->is_leaf);
        for (auto& [c, link] : loc_curr->links) {
            auto start = link->start;
            auto end = link->end;
            auto* loc_child = link->node.get();
            auto len = end - start + 1;
            if (reverse) {
                len = start - end + 1;
            }
            node->links[c] = std::make_unique<Blind_tree_link>(len, dfs(loc_child));
        }
        return node.release();
    };

    this->root = std::unique_ptr<Node>(dfs(location_tree.root.get()));
}

// adjust match to match the blind tree node and links
// std::vector<size_t> Blind_tree::match(const std::string& pattern) {
std::vector<size_t> Blind_tree::match(std::string_view pattern, size_t factor_size, size_t offset) {
    if (!reverse) {
        auto curr = root.get();
        size_t i = 0;
        while (i < pattern.size() && curr->links.size()) {
            auto it = curr->links.find(pattern[i]);
            if (it == curr->links.end()) {
                break;
            }
            auto len = it->second->len;
            auto& node = it->second->node;
            // if (i + len - 1 < pattern.size()) {                
            if (i + len >= factor_size) {
                // get all locations
                std::vector<size_t> result;
                get_locations(node.get(), result, offset);

                // // pick one location and verify
                // auto location = result[0];
                // auto pattern_in_text = std::string_view(text).substr(location, pattern.size());
                // if (pattern != pattern_in_text) {
                //     return {};
                // }

                return result;
            }
            curr = node.get();
            // i = i + len - 1;
            i = i + len;
        }
        return {};
    }

    // reverse
    auto curr = root.get();
    long long i = pattern.size() - 1; // long long to avoid underflow
    while (i >= 0 && curr->links.size()) {
        auto it = curr->links.find(pattern[i]);
        if (it == curr->links.end()) {
            break;
        }
        auto len = it->second->len;
        auto& node = it->second->node;
        // if (i - len + 1 >= 0) {
        // if ((i - len) <= (pattern.size() - 1 - factor_size)) { // I don't know why this format won't work
        // if ((pattern.size() - i + len - 1) >= factor_size) {
        if (i < len) {
            // get all locations
            std::vector<size_t> result;
            get_locations(node.get(), result, offset);

            // // pick one location and verify
            // auto location = result[0];
            // auto pattern_in_text = std::string_view(text).substr(location + factor_size - pattern.size(), pattern.size());
            // if (pattern != pattern_in_text) {
            //     return {};
            // }

            return result;
        }
        curr = node.get();
        // i = i - len + 1;
        i = i - len;
    }
    return {};
}

// get number of nodes
size_t Blind_tree::get_num_nodes() {
    size_t result = 0;
    std::stack<Node*> s;
    s.push(root.get());
    while (s.size()) {
        auto curr = s.top();
        s.pop();
        result++;
        for (auto& [c, link] : curr->links) {
            s.push(link->node.get());
        }
    }
    return result;
}

void Blind_tree::serialize(std::ofstream& out) {
    std::vector<bool> LOUDS_tmp, is_leaf_tmp;
    std::vector<size_t> label_tmp, factor_idx_tmp, len_tmp;

    std::queue<Node*> que;
    que.push(root.get());
    while (!que.empty()) {
        Node* cur = que.front();
        que.pop();
        LOUDS_tmp.push_back(1);
        
        for (const auto &p : cur->links) {
            label_tmp.emplace_back(p.first);
            len_tmp.emplace_back(p.second->len);
            que.push(p.second->node.get());
            LOUDS_tmp.push_back(0);
        }
        
        is_leaf_tmp.emplace_back(cur->is_leaf);
        factor_idx_tmp.emplace_back(cur->factor_index);
    }

    sdsl::bit_vector LOUDS, is_leaf;
    sdsl::int_vector<> label, factor_idx, len;

    LOUDS.resize(LOUDS_tmp.size());
    for (size_t i = 0, e = LOUDS_tmp.size(); i < e; i++) { LOUDS[i] = LOUDS_tmp[i]; }
    is_leaf.resize(is_leaf_tmp.size());
    for (size_t i = 0, e = is_leaf_tmp.size(); i < e; i++) { is_leaf[i] = is_leaf_tmp[i]; }
    label.resize(label_tmp.size());
    for (size_t i = 0, e = label_tmp.size(); i < e; i++) { label[i] = label_tmp[i]; }
    factor_idx.resize(factor_idx_tmp.size());
    for (size_t i = 0, e = factor_idx_tmp.size(); i < e; i++) { factor_idx[i] = factor_idx_tmp[i]; }
    len.resize(len_tmp.size());
    for (size_t i = 0, e = len_tmp.size(); i < e; i++) { len[i] = len_tmp[i]; }

    sdsl::util::bit_compress(LOUDS);
    sdsl::util::bit_compress(is_leaf);
    sdsl::util::bit_compress(label);
    sdsl::util::bit_compress(factor_idx);
    sdsl::util::bit_compress(len);

    LOUDS.serialize(out);
    is_leaf.serialize(out);
    label.serialize(out);
    factor_idx.serialize(out);
    len.serialize(out);

}

void Blind_tree::load(std::ifstream& in) {
    sdsl::bit_vector LOUDS, is_leaf;
    sdsl::int_vector<> label, factor_idx, len;
    size_t LOUDS_idx = 0, is_leaf_idx = 0, label_idx = 0, factor_idx_idx = 0, len_idx = 0;
    
    LOUDS.load(in);
    is_leaf.load(in);
    label.load(in);
    factor_idx.load(in);
    len.load(in);

    auto read_node_info = [&](Node* cur) {
        if (cur == nullptr) { throw std::invalid_argument("false_0"); }
        if (is_leaf_idx >= is_leaf.size()) { throw std::invalid_argument("false_1"); }
        cur->is_leaf = is_leaf[is_leaf_idx++];
        if (factor_idx_idx >= factor_idx.size()) { throw std::invalid_argument("false_2"); }
        cur->factor_index = factor_idx[factor_idx_idx++];
    };

    std::queue<Node*> que;
    root = std::make_unique<Node>();
    Node* cur = root.get();
    read_node_info(cur);
    for (size_t i = 1, e = LOUDS.size(); i < e; i++) {
        if (LOUDS[i] == 0) {
            std::unique_ptr new_node = std::make_unique<Node>();
            que.push(new_node.get());
            cur->links[(char) label[label_idx++]] = std::make_unique<Link>(len[len_idx++], new_node.release());
        } else {
            cur = que.front();
            que.pop();
            read_node_info(cur);
        }
    }

}

void Blind_tree::print(Node* node, bool first_call, bool print_out_path) {
    if (first_call) {
        node = root.get();
    }
    if (!node) {
        return;
    }
    std::cout << "node:" << reinterpret_cast<size_t>(node) << std::endl;
    std::cout << "factor_index: " << node->factor_index << std::endl;
    std::cout << "is_leaf: " << node->is_leaf << std::endl;
    std::cout << "links: " << std::endl;
    for (auto& [c, link] : node->links) {
        std::cout << "\tc: " << c << std::endl;
        std::cout << "\tlen: " << link->len << std::endl;
        std::cout << "\tnode: " << reinterpret_cast<size_t>(link->node.get()) << std::endl;
    }
    std::cout << std::endl;
    for (auto& [c, link] : node->links) {
        print(link->node.get(), false, print_out_path);
    }
}

// move the get_locations function out of the function
void Blind_tree::get_locations(Node* node, std::vector<size_t>& result, size_t offset) {
    if (!node) {
        return;
    }
    if (node->is_leaf) {
        result.push_back(node->factor_index - offset);
    }
    for (auto& [c, link] : node->links) {
        get_locations(link->node.get(), result,  offset);
    }
} 

size_t Blind_tree::pick_one_location(Node* node) {
    if (!node->links.size() || node->is_leaf) {
        return node->factor_index;
    }

    // pick one and iterate
    for (auto& [c, link] : node->links) {
        if (pick_one_location(link->node.get()) != -1) {
            return pick_one_location(link->node.get());
        }
    }

    return -1;
}