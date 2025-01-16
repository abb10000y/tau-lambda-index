// Compact trie structure for location matching
// TODO combine the forward and reverse version

#pragma once
#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <map>
#include <functional> // TODO remove
#include "Location_tree_node.h"
#include "Location_tree_link.h"

class Location_tree {
public:
    using Node = Location_tree_node;
    using Link = Location_tree_link;

    // construct a location tree from the text
    // reverse: whether to construct a reverse location tree
    Location_tree(const std::string& text, bool reverse = false);
    ~Location_tree();

    void insert_factors(const std::vector<std::pair<size_t, size_t>>& min_factors);
    std::vector<size_t> match(const std::string& pattern);

    // debug, print info of the tree for each node
    void print(Node* node = nullptr, bool first_call = true, bool print_out_path = true); 
    size_t get_num_nodes(); // TODO remove

    void serialize(std::string& filename);
    void load(std::string& filename);

    std::unique_ptr<Node> root; // Blind tree needs to access this, TODO better way to do this?

private:
    // std::unique_ptr<Node> root;
    const std::string& text;
    bool reverse;
    bool debug = false; // TODO remove

    // insert one factor into the tree
    void foward_insert_factor(const std::pair<size_t, size_t>& min_factor);
    void reverse_insert_factor(const std::pair<size_t, size_t>& min_factor);

    // util for foward and reverse matching
    std::vector<size_t> foward_match(const std::string& pattern);
    std::vector<size_t> reverse_match(const std::string& pattern);

    // recursively get all the locations of the subtree
    // TODO non dfs version
    void get_locations(Node* node, std::vector<size_t>& result); 
};

Location_tree::Location_tree(const std::string& text, bool reverse) 
    : text(text), root(std::make_unique<Node>()), reverse(reverse) {}

Location_tree::~Location_tree() {}

void Location_tree::insert_factors(const std::vector<std::pair<size_t, size_t>>& min_factors) {
    if (reverse) {
        for (const auto& min_factor : min_factors) {
            reverse_insert_factor(min_factor);
        }
    } else {
        // reversely insert the factors, to prevent the case that a is prefix of b and b is inserted first
        // TODO better way to do this?
        for (auto itr = min_factors.rbegin(); itr != min_factors.rend(); ++itr) {
            foward_insert_factor(*itr);
        }
    }
    if (debug) { // TODO remove
        print();
    }
}

void Location_tree::foward_insert_factor(const std::pair<size_t, size_t>& min_factor) {
    size_t start = min_factor.first;
    size_t end = text.size() - 1;
    size_t factor_index = min_factor.first;

    // find the node to insert
    Node* node = root.get();
    while (start <= end) {
        char c = text[start];
        if (node->links.find(c) == node->links.end()) {
            node->links[c] = std::make_unique<Link>(start, end, factor_index);
            return;
        }

        // pick one link and traverse
        Link* link = node->links[c].get();
        size_t i = link->start;
        size_t j = link->end;
        while (i <= j && start <= end) {
            if (text[i] != text[start]) {
                break;
            }
            ++i;
            ++start;
        }

        // continue if the link is exhausted
        if (i > j) {
            node = link->node.get();
            continue;
        }

        // split the link if the link is not exhausted
        auto split_point = std::make_unique<Node>();
        split_point->links[text[i]] = std::make_unique<Link>(i, j, std::move(link->node)); // TODO move twice?
        split_point->links[text[start]] = std::make_unique<Link>(start, end, factor_index);
        link->end = i - 1;
        link->node = std::move(split_point);
        return;
    }
}

// TODO more clean way to implement this
// TODO one case might cause problem: when a is prefix of b and b is inserted first
void Location_tree::reverse_insert_factor(const std::pair<size_t, size_t>& min_factor) {
    long long start = min_factor.second; // TODO to prevent underflow, better way to do this?
    long long end = 0;
    size_t factor_index = min_factor.first;

    // TODO remove
    if (debug) {
        std::cout << "reverse insert factor: " << text.substr(end , start - end + 1);
        std::cout << " left: " << min_factor.first << " right: " << min_factor.second << std::endl;
    }

    // find the node to insert
    Node* node = root.get();
    while (start >= end) {
        char c = text[start];
        if (node->links.find(c) == node->links.end()) {
            node->links[c] = std::make_unique<Link>(start, end, factor_index);
            if (debug) std::cout << "create node: " << c << std::endl;
            return;
        }

        // pick one link and traverse
        Link* link = node->links[c].get();
        long long i = link->start; // TODO to prevent underflow, better way to do this?
        long long j = link->end;
        while (i >= j && start >= end) {
            if (text[i] != text[start]) {
                break;
            }
            if (debug) std::cout << "i: " << i << " start: " << start << std::endl;
            --i;
            --start;
        }
        
        if (debug) {
            std::cout << "itr done: i = " << i << " j = " << j << " start = " << start << " end = " << end << std::endl;
            std::cout << ((i < static_cast<long long> (j)) ? "i < j" : "i >= j") << std::endl;
        }

        // continue if the link is exhausted
        if (i < j) {
            node = link->node.get();
            if (debug) std::cout << "i < j: " << i << " " << j << std::endl;
            continue;
        }

        if (debug) std::cout << "split" << std::endl;
        // split the link if the link is not exhausted
        auto split_point = std::make_unique<Node>();
        split_point->links[text[i]] = std::make_unique<Link>(i, j, std::move(link->node)); // TODO move twice?
        split_point->links[text[start]] = std::make_unique<Link>(start, end, factor_index);
        link->end = i + 1;
        link->node = std::move(split_point);
        return;
    }
}

// find the first node that matches the pattern
std::vector<size_t> Location_tree::match(const std::string& pattern) {
    if (reverse) {
        return reverse_match(pattern);
    } else {
        return foward_match(pattern);
    }
}

std::vector<size_t> Location_tree::foward_match(const std::string& pattern) {
    Node* node = root.get();
    size_t i = 0;
    while (i < pattern.size()) {
        char c = pattern[i];
        if (node->links.find(c) == node->links.end()) {
            return {};
        }

        Link* link = node->links[c].get();
        size_t j = link->start;
        while (i < pattern.size() && j <= link->end) {
            if (pattern[i] != text[j]) {
                return {};
            }
            ++i;
            ++j;
        }

        if (i == pattern.size()) {
            // found a match
            std::vector<size_t> result;
            get_locations(link->node.get(), result);
            return result;
        }

        node = link->node.get();
    }

    return {};
}

// Current implementation is not correct TODO
std::vector<size_t> Location_tree::reverse_match(const std::string& pattern) {
    Node* node = root.get();
    long long i = pattern.size() - 1; // TODO to prevent underflow, better way to do this?

    if (debug) {
        std::cout << "\n\nreverse match: " << pattern << std::endl;
    }

    while (i >= 0) {
        char c = pattern[i];
        if (node->links.find(c) == node->links.end()) {
            if (debug) std::cout << "not found: " << c << std::endl;
            return {};
        }

        Link* link = node->links[c].get();
        long long j = link->start; // TODO to prevent underflow, better way to do this?

        if (debug) std::cout << "link: " << c << " start: " << link->start << " end: " << link->end << std::endl;

        while (i >= 0 && j >= static_cast<long long> (link->end)) { // TODO cast to prevent underflow, remove?
            if (debug) std::cout << "i: " << i << " j: " << j << " pattern[i]: " << pattern[i] << " text[j]: " << text[j] << std::endl;
            if (pattern[i] != text[j]) {
                if (debug) std::cout << "not match" << std::endl;
                return {};
            }
            --i;
            --j;
        }

        if (debug) {
            std::cout << "itr done: i = " << i << " j = " << j << std::endl;
            std::cout << ((i < static_cast<long long> (j)) ? "i < j" : "i >= j") << std::endl;
        }

        if (i < 0) { // TODO index bound?
            // found a match
            if (debug) std::cout << "found a match" << std::endl;
            std::vector<size_t> result;
            get_locations(link->node.get(), result);
            return result;
        }

        if (debug) std::cout << "continue" << std::endl;
        node = link->node.get();
    }

    if (debug) std::cout << "return" << std::endl;
    return {};
}

// utility function to get all the locations of the subtree
void Location_tree::get_locations(Node* node, std::vector<size_t>& result) {
    if (node->is_leaf) {
        result.push_back(node->factor_index);
    }

    for (const auto& [c, link] : node->links) {
        get_locations(link->node.get(), result);
    }
}

// debug, print info of the tree for each node
void Location_tree::print(Node* node, bool first_call, bool print_out_path) {
    if (first_call) {
        node = root.get();
        std::cout << "The text is: " << text << "\n\n";
    }
    std::cout << "node: " << reinterpret_cast<size_t>(node) << std::endl;
    std::cout << "factor_index: " << node->factor_index << std::endl;
    std::cout << "is_leaf: " << node->is_leaf << std::endl;
    // info for each link
    for (const auto& [c, link] : node->links) {
        std::cout << "\tlink: " << std::endl;
        // std::cout << "\tstart: " << link->start << std::endl;
        // std::cout << "\tend: " << link->end << std::endl;
        
        // debug, first char and len
        std::cout << "\tfirst char: " << text[link->start] << std::endl;
        std::cout << "\tlen: " << link->end - link->start + 1 << std::endl;

        // print path labels
        std::cout << "\tpath: ";
        if (print_out_path) {
            if (reverse) {
                for (long long i = link->start; i >= (long long)(link->end); --i) {
                    std::cout << text[i];
                }
            } else {
                for (size_t i = link->start; i <= link->end; ++i) {
                    std::cout << text[i];
                }
            }
        }

        std::cout << std::endl;
        std::cout << "\tnode: " << reinterpret_cast<size_t>(link->node.get()) << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;

    for (const auto& [c, link] : node->links) {
        print(link->node.get(), false);
    }
}

size_t Location_tree::get_num_nodes() {
    std::function<size_t(Node*)> dfs = [&](Node* node) {
        if (node == nullptr) {
            return static_cast<size_t>(0);
        }

        size_t result = 1;
        for (const auto& [c, link] : node->links) {
            result += dfs(link->node.get());
        }
        return result;
    };

    return dfs(root.get());
}

void Location_tree::serialize(std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    using index_type = uint32_t;

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

        out.write(reinterpret_cast<char*>(&id), sizeof(id));
        // out.write(reinterpret_cast<char*>(&is_leaf), sizeof(is_leaf));
        // out.write(reinterpret_cast<char*>(&children_num), sizeof(children_num));
        out.write(reinterpret_cast<char*>(&is_leaf_children_num), sizeof(is_leaf_children_num));
        if (is_leaf) {
            out.write(reinterpret_cast<char*>(&factor_index), sizeof(factor_index));
        } 
            
        // use char because the num of child is always fewer than 128
        for (const auto& p : node->links) {
            char c = p.first;
            Link *link = p.second.get();
            index_type start = link->start;
            index_type end = link->end;
            index_type child_id = node_id[link->node.get()];
            out.write(reinterpret_cast<char*>(&c), sizeof(c));
            out.write(reinterpret_cast<char*>(&start), sizeof(start));
            out.write(reinterpret_cast<char*>(&end), sizeof(end));
            out.write(reinterpret_cast<char*>(&child_id), sizeof(child_id));
        }
    }
}

void Location_tree::load(std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    using index_type = uint32_t;

    // id_node mappping, and id_link mapping, for first pass
    std::unordered_map<index_type, Node*> id_node;
    std::unordered_map<index_type, std::vector<std::tuple<char, index_type, index_type, index_type>>> id_links; // tuple<c, start, end, child_id>

    index_type node_num;
    in.read(reinterpret_cast<char*>(&node_num), sizeof(node_num));

    // first pass
    for (index_type i = 0; i < node_num; ++i) {
        index_type id;
        bool is_leaf;
        uint8_t children_num;
        index_type factor_index;

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

        for (int j = 0; j < children_num; ++j) {
            char c;
            index_type start;
            index_type end;
            index_type child_id;
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            in.read(reinterpret_cast<char*>(&start), sizeof(start));
            in.read(reinterpret_cast<char*>(&end), sizeof(end));
            in.read(reinterpret_cast<char*>(&child_id), sizeof(child_id));
            if (id_links.find(id) == id_links.end()) {
                id_links[id] = std::vector<std::tuple<char, index_type, index_type, index_type>>();
            }
            id_links[id].push_back(std::make_tuple(c, start, end, child_id));
        }
    }
    
    // rewrite the buggy code above again
    for (auto& p : id_links) {
        index_type id = p.first;
        for (auto& link : p.second) {
            char c;
            index_type start;
            index_type end;
            index_type child_id;
            std::tie(c, start, end, child_id) = link;
            id_node[id]->links[c] = std::make_unique<Link>(start, end, std::unique_ptr<Node>(id_node[child_id]));
        }
    }

    root = std::unique_ptr<Node>(id_node[0]);
}