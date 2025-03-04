// Implementation of the Ac_automata class.
// For min-factor matching.

#pragma once

#include <string>
#include <vector>
#include <map>
#include <queue>
#include <tuple>
#include <stack>
#include <fstream>
#include <memory>
#include <sdsl/bit_vectors.hpp>


class Ac_automata {
public:
    Ac_automata() = default; // Add default constructor to make Index class work, TODO better solution?
    Ac_automata(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors);
    ~Ac_automata();

    // move assginment operator, to make Index class work, TODO better solution?
    Ac_automata& operator=(Ac_automata&& other) {
        // root = std::move(other.root);
        node_vector = std::move(other.node_vector);
        return *this;
    }

    // return the positions of the matched min-factors
    // return type: vector of tuples, each tuple contains the start_pos_in_text, start_pos_in_pattern, length
    // TODO better return type
    std::vector<std::tuple<size_t, size_t, size_t>> match(const std::string& pattern, bool match_all); 

    // return <pos_in_pat, length>
    std::pair<size_t, size_t> match_pos_in_pattern(const std::string& pattern);

    // debugging info TODO remove
    size_t get_ac_size() {
        return node_vector.size();
    }

    void serialize(std::ofstream &out);
    void load(std::ifstream &in);

private:
    // define the node of the trie
    struct Node {
        std::map<char, size_t> children;
        // Node* fail; // pointer to the failure node
        size_t fail = 0; // node_vector.size() as nullptr

        bool is_min_factor = false;
        std::vector<size_t> positions; // positions of the matched min-factors in the text
        size_t height = 0; // height of the node in the trie, i.e. the length of the string represented by the node
    };

    // std::unique_ptr<Node> root;
    size_t root = 0;
    std::vector<std::unique_ptr<Node>> node_vector;
    void construct_trie(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors);
    void construct_failure_links();
};

Ac_automata::Ac_automata(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors) {
    std::cout << "Constructing AC_automaton" << std::endl;
    construct_trie(text, min_factors);
    construct_failure_links();
}

Ac_automata::~Ac_automata() {
}

void Ac_automata::construct_trie(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors) {
    node_vector.push_back(std::make_unique<Node>());
    for (auto p : min_factors) {
        std::string_view factor = std::string_view(text).substr(p.first, p.second - p.first + 1);
        size_t cur = root;
        for (char c : factor) {
            if (node_vector[cur]->children.find(c) == node_vector[cur]->children.end()) {
                node_vector[cur]->children[c] = node_vector.size();
                node_vector.emplace_back(std::make_unique<Node>());
            }
            cur = node_vector[cur]->children[c];
        }
        node_vector[cur]->is_min_factor = true;
        node_vector[cur]->positions.emplace_back(p.first);
        node_vector[cur]->height = factor.size();
    }
}

void Ac_automata::construct_failure_links() {
    size_t null_link = node_vector.size();

    node_vector[0]->fail = null_link;
    std::queue<size_t> q;
    q.push(root);
    while (!q.empty()) {
        size_t cur = q.front();
        q.pop();
        for (const auto& p : node_vector[cur]->children) {
            char c = p.first;
            size_t child = p.second;
            size_t fail_node = node_vector[cur]->fail;
            while (fail_node != null_link && node_vector[fail_node]->children.find(c) == node_vector[fail_node]->children.end()) {
                fail_node = node_vector[fail_node]->fail;
            }
            if (fail_node == null_link) {
                node_vector[child]->fail = root;
            } else {
                node_vector[child]->fail = node_vector[fail_node]->children[c];
            }
            q.push(child);
        }
    }
}

std::vector<std::tuple<size_t, size_t, size_t>> Ac_automata::match(const std::string& pattern, bool match_all = false) {
    std::vector<std::tuple<size_t, size_t, size_t>> matched_min_factors;
    size_t cur = root, null_link = node_vector.size();
    for (size_t i = 0; i < pattern.size(); i++) {
        char c = pattern[i];
        
        // find the next node
        while (cur != null_link && node_vector[cur]->children.find(c) == node_vector[cur]->children.end()) {
            cur = node_vector[cur]->fail;
        }

        if (cur == null_link) {
            cur = root;
            continue;
        }

        cur = node_vector[cur]->children[c];
        
        if (node_vector[cur]->is_min_factor) {
            for (size_t pos : node_vector[cur]->positions) {
                matched_min_factors.emplace_back(pos, i - node_vector[cur]->height + 1, node_vector[cur]->height);
            }
            if (!match_all) {
                break;
            }
        }
    }
    return matched_min_factors;   
}

std::pair<size_t, size_t> Ac_automata::match_pos_in_pattern(const std::string& pattern) {
    std::vector<size_t> matched_min_factors_pos_in_text;
    size_t cur = root, null_link = node_vector.size();
    for (size_t i = 0; i < pattern.size(); i++) {
        char c = pattern[i];
        
        // find the next node
        while (cur != null_link && node_vector[cur]->children.find(c) == node_vector[cur]->children.end()) {
            cur = node_vector[cur]->fail;
        }

        if (cur == null_link) {
            cur = root;
            continue;
        }

        cur = node_vector[cur]->children[c];
        
        if (node_vector[cur]->is_min_factor) {
            return std::make_pair(i - node_vector[cur]->height + 1, node_vector[cur]->height);
        }
    }
    return std::make_pair(-1, -1);   
}

// output in binary format
void Ac_automata::serialize(std::ofstream &out) {
    std::unordered_map<size_t, size_t> nodeVectorIdx_to_levelOrder;
    std::vector<bool> LOUDS_tmp, is_min_factor_tmp;
    std::vector<size_t> label_tmp, position_tmp, position_cnt_tmp, height_tmp, failure_tmp;
    size_t null_link = node_vector.size(), level_idx = 0;
    nodeVectorIdx_to_levelOrder[null_link] = null_link;
    
    // transform into LOUDS representation
    std::queue<size_t> que;
    que.push(root);
    while (!que.empty()) {
        size_t cur = que.front();
        que.pop();
        LOUDS_tmp.emplace_back(1);
        nodeVectorIdx_to_levelOrder[cur] = level_idx++;

        for (const auto &p : node_vector[cur]->children) {
            LOUDS_tmp.emplace_back(0);
            label_tmp.emplace_back(p.first);
            que.push(p.second);
        }

        failure_tmp.emplace_back(nodeVectorIdx_to_levelOrder[node_vector[cur]->fail]);

        is_min_factor_tmp.emplace_back(node_vector[cur]->is_min_factor);
        
        position_cnt_tmp.emplace_back(node_vector[cur]->positions.size());
        for (const auto p : node_vector[cur]->positions) { position_tmp.emplace_back(p); }
        
        height_tmp.emplace_back(node_vector[cur]->height);
    }

    // transform into sdsl data structure
    sdsl::bit_vector LOUDS, is_min_factor;
    sdsl::int_vector<> label, position, position_cnt, height, failure;
    LOUDS.resize(LOUDS_tmp.size());
    for (size_t i = 0, e = LOUDS_tmp.size(); i < e; i++) { LOUDS[i] = LOUDS_tmp[i]; }
    is_min_factor.resize(is_min_factor_tmp.size());
    for (size_t i = 0, e = is_min_factor_tmp.size(); i < e; i++) { is_min_factor[i] = is_min_factor_tmp[i]; }
    label.resize(label_tmp.size());
    for (size_t i = 0, e = label_tmp.size(); i < e; i++) { label[i] = label_tmp[i]; }
    position.resize(position_tmp.size());
    for (size_t i = 0, e = position_tmp.size(); i < e; i++) { position[i] = position_tmp[i]; }
    position_cnt.resize(position_cnt_tmp.size());
    for (size_t i = 0, e = position_cnt_tmp.size(); i < e; i++) { position_cnt[i] = position_cnt_tmp[i]; }
    height.resize(height_tmp.size());
    for (size_t i = 0, e = height_tmp.size(); i < e; i++) { height[i] = height_tmp[i]; }
    failure.resize(failure_tmp.size());
    for (size_t i = 0, e = failure_tmp.size(); i < e; i++) { failure[i] = failure_tmp[i]; }
    
    sdsl::util::bit_compress(LOUDS);
    sdsl::util::bit_compress(is_min_factor);
    sdsl::util::bit_compress(label);
    sdsl::util::bit_compress(position);
    sdsl::util::bit_compress(position_cnt);
    sdsl::util::bit_compress(height);
    sdsl::util::bit_compress(failure);


    // serializing
    sdsl::write_member(null_link, out);
    LOUDS.serialize(out);
    is_min_factor.serialize(out);
    label.serialize(out);
    position.serialize(out);
    position_cnt.serialize(out);
    height.serialize(out);
    failure.serialize(out);
}

// read in binary format
void Ac_automata::load(std::ifstream &in) {
    sdsl::bit_vector LOUDS, is_min_factor;
    sdsl::int_vector<> label, position, position_cnt, height, failure;
    size_t node_idx = 1, is_min_factor_idx = 0, label_idx = 0, position_idx = 0, position_cnt_idx = 0, height_idx = 0, failure_idx = 0, null_link;

    sdsl::read_member(null_link, in);
    LOUDS.load(in);
    is_min_factor.load(in);
    label.load(in);
    position.load(in);
    position_cnt.load(in);
    height.load(in);
    failure.load(in);

    node_vector.resize(null_link);

    std::cout << "(TBD) ac_size: " << null_link << std::endl;
    std::cout << "(TBD) LOUDS.size: " << LOUDS.size() << std::endl;
    
    auto read_node_info = [&](size_t cur, size_t idx) {
        // if (failure[failure_idx] != null_link && failure[failure_idx] >= idx) { throw std::invalid_argument("a"); }
        node_vector[cur]->fail = failure[failure_idx++];

        node_vector[cur]->is_min_factor = is_min_factor[is_min_factor_idx++];

        for (size_t i = 0, e = position_cnt[position_cnt_idx++]; i < e; i++) { node_vector[cur]->positions.emplace_back(position[position_idx++]); }
        
        node_vector[cur]->height = height[height_idx++];
    };
    
    node_vector[root] = std::make_unique<Node>();
    size_t cur = root, idx = 1;
    read_node_info(cur, idx);
    for (size_t i = 1, e = LOUDS.size(); i < e; i++) {
        if (LOUDS[i] == 0) {
            node_vector[cur]->children[(char) label[label_idx++]] = idx;
            node_vector[idx++] = std::make_unique<Node>();
        } else {
            cur++;
            read_node_info(cur, idx);
        }
    }
}