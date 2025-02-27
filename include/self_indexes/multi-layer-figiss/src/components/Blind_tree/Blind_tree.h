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
#include <functional> // TODO remove
#include "../Location_tree/Location_tree.h"
#include "Blind_tree_node.h"
#include "Blind_tree_link.h"

class Blind_tree {
public:
    using Node = Blind_tree_node;

    // construct a location tree from the text
    // reverse: whether to construct a reverse location tree
    Blind_tree(std::string& text, bool reverse = false);
    ~Blind_tree();

    void insert_factors(const std::vector<std::pair<size_t, size_t>>& min_factors);
    // std::vector<size_t> match(const std::string& pattern);
    std::vector<size_t> match(std::string_view pattern, size_t factor_size);

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
    void get_locations(Node* node, std::vector<size_t>& result); 

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
std::vector<size_t> Blind_tree::match(std::string_view pattern, size_t factor_size) {
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
            if (i + len + 1 >= factor_size) {
                // get all locations
                std::vector<size_t> result;
                get_locations(node.get(), result);

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
        if ((pattern.size() - i + len) >= factor_size) {
            // get all locations
            std::vector<size_t> result;
            get_locations(node.get(), result);

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
    using index_type = uint32_t;

    out.write(reinterpret_cast<char*>(&reverse), sizeof(reverse));

    std::unordered_map<Node*, index_type> node_id;
    // dfs to assign id to each node
    std::stack<Node*> s;
    s.push(root.get());
    index_type id = 0;
    while (!s.empty()) {
        Node* curr_node = s.top();
        s.pop();
        node_id[curr_node] = id++;
        for (const auto& p : curr_node->links) {
            s.push(p.second->node.get());
        }
    }

    // serialize
    index_type node_num = node_id.size();
    out.write(reinterpret_cast<char*>(&node_num), sizeof(node_num));
    for (const auto& p: node_id) {
        Node* node = p.first;

        index_type id = p.second;
        bool is_leaf = node->is_leaf; // optimize
        index_type factor_index = node->factor_index;
        uint8_t children_num = node->links.size(); // optimize

        // combine is_leaf and children_num into one byte
        uint8_t is_leaf_children_num = (is_leaf << 7) | children_num;

        // write primitive vars
        out.write(reinterpret_cast<char*>(&id), sizeof(id));
        // out.write(reinterpret_cast<char*>(&is_leaf), sizeof(is_leaf));
        // out.write(reinterpret_cast<char*>(&children_num), sizeof(children_num));
        out.write(reinterpret_cast<char*>(&is_leaf_children_num), sizeof(is_leaf_children_num));
        if (is_leaf) {
            out.write(reinterpret_cast<char*>(&factor_index), sizeof(factor_index));
        }

        
            
        for (const auto& p : node->links) {
            Blind_tree_link *link = p.second.get();

            char c = p.first; // optimize
            u_int16_t len = link->len; // bound by lambda
            index_type child_id = node_id[link->node.get()];

            out.write(reinterpret_cast<char*>(&c), sizeof(c));
            out.write(reinterpret_cast<char*>(&len), sizeof(len));
            out.write(reinterpret_cast<char*>(&child_id), sizeof(child_id));
        }
    }
}

void Blind_tree::load(std::ifstream& in) {
    using index_type = uint32_t;

    in.read(reinterpret_cast<char*>(&reverse), sizeof(reverse));

    // id_node, id_link mapping
    std::unordered_map<index_type, Node*> id_node;
    std::unordered_map<index_type, std::vector<std::tuple<char, u_int16_t, index_type>>> id_link;

    // deserialize
    index_type node_num;
    in.read(reinterpret_cast<char*>(&node_num), sizeof(node_num));

    // first pass, create nodes
    for (index_type i = 0; i < node_num; i++) {

        index_type id;
        bool is_leaf; // optimize
        index_type factor_index;
        uint8_t children_num; // optimize

        // combined vars
        uint8_t is_leaf_children_num;

        in.read(reinterpret_cast<char*>(&id), sizeof(id));
        // in.read(reinterpret_cast<char*>(&is_leaf), sizeof(is_leaf));
        // in.read(reinterpret_cast<char*>(&children_num), sizeof(children_num));
        in.read(reinterpret_cast<char*>(&is_leaf_children_num), sizeof(is_leaf_children_num));
        is_leaf = is_leaf_children_num >> 7;
        children_num = is_leaf_children_num & 0x7f;

        id_node[id] = new Node();
        id_node[id]->is_leaf = is_leaf;
        if (is_leaf) {
            in.read(reinterpret_cast<char*>(&factor_index), sizeof(factor_index));
            id_node[id]->factor_index = factor_index;
        }

        for (int j = 0; j < children_num; j++) {

            char c; // optimize
            u_int16_t len; // bound by lambda
            index_type child_id;

            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            in.read(reinterpret_cast<char*>(&len), sizeof(len));
            in.read(reinterpret_cast<char*>(&child_id), sizeof(child_id));
            if (id_link.find(id) == id_link.end()) {
                id_link[id] = {};
            }
            id_link[id].emplace_back(c, len, child_id);
        }
    }

    // second pass, create links
    for (const auto& p : id_link) {
        index_type id = p.first;
        for (const auto& [c, len, child_id] : p.second) {
            id_node[id]->links[c] = std::make_unique<Blind_tree_link>(len, id_node[child_id]);
        }
    }

    this->root = std::unique_ptr<Node>(id_node[0]);
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
void Blind_tree::get_locations(Node* node, std::vector<size_t>& result) {
    if (!node) {
        return;
    }
    if (node->is_leaf) {
        result.push_back(node->factor_index);
    }
    for (auto& [c, link] : node->links) {
        get_locations(link->node.get(), result);
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