// Implementation of the Ac_automata class.
// For min-factor matching.

#pragma once

#include <string>
#include <vector>
#include <map>
#include <queue>
#include <tuple>
#include <memory>

class Ac_automata {
public:
    Ac_automata() = default; // Add default constructor to make Index class work, TODO better solution?
    Ac_automata(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors);
    ~Ac_automata();

    // move assginment operator, to make Index class work, TODO better solution?
    Ac_automata& operator=(Ac_automata&& other) {
        root = std::move(other.root);
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
        return ac_size;
    }

    void serialize(std::string& filename);
    void load(std::string& filename);

private:
    // define the node of the trie
    struct Node {
        std::map<char, std::unique_ptr<Node>> children;
        Node* fail; // pointer to the failure node

        bool is_min_factor;
        std::vector<size_t> positions; // positions of the matched min-factors in the text
        size_t height; // height of the node in the trie, i.e. the length of the string represented by the node
    };

    std::unique_ptr<Node> root;
    void construct_trie(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors);
    void construct_failure_links();

    // debugging info TODO remove
    size_t ac_size;
};

Ac_automata::Ac_automata(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors) {
    construct_trie(text, min_factors);
    construct_failure_links();
}

Ac_automata::~Ac_automata() {
}

void Ac_automata::construct_trie(const std::string& text, const std::vector<std::pair<size_t, size_t>>& min_factors) {
    ac_size = 0;
    root = std::make_unique<Node>();
    for (auto p : min_factors) {
        std::string_view factor = std::string_view(text).substr(p.first, p.second - p.first + 1);
        Node* curr_node = root.get();
        for (char c : factor) {
            if (curr_node->children.find(c) == curr_node->children.end()) {
                curr_node->children[c] = std::make_unique<Node>();
                ++ac_size;
            }
            curr_node = curr_node->children[c].get();
        }
        curr_node->is_min_factor = true;
        curr_node->positions.emplace_back(p.first);
        curr_node->height = factor.size();
    }
}

void Ac_automata::construct_failure_links() {
    root->fail = nullptr;
    std::queue<Node*> q;
    q.push(root.get());
    while (!q.empty()) {
        Node* curr_node = q.front();
        q.pop();
        for (const auto& p : curr_node->children) {
            char c = p.first;
            Node* child = p.second.get();
            Node* fail_node = curr_node->fail;
            while (fail_node != nullptr && fail_node->children.find(c) == fail_node->children.end()) {
                fail_node = fail_node->fail;
            }
            if (fail_node == nullptr) {
                child->fail = root.get();
            } else {
                child->fail = fail_node->children[c].get();
            }
            q.push(child);
        }
    }
}

std::vector<std::tuple<size_t, size_t, size_t>> Ac_automata::match(const std::string& pattern, bool match_all = false) {
    std::vector<std::tuple<size_t, size_t, size_t>> matched_min_factors;
    Node* curr_node = root.get();
    for (size_t i = 0; i < pattern.size(); i++) {
        char c = pattern[i];
        
        // find the next node
        while (curr_node != nullptr && curr_node->children.find(c) == curr_node->children.end()) {
            curr_node = curr_node->fail;
        }

        if (curr_node == nullptr) {
            curr_node = root.get();
            continue;
        }

        curr_node = curr_node->children[c].get();
        
        if (curr_node->is_min_factor) {
            for (size_t pos : curr_node->positions) {
                matched_min_factors.emplace_back(pos, i - curr_node->height + 1, curr_node->height);
            }
            if (!match_all) {
                break;
            }
        }
    }
    return matched_min_factors;
}

std::pair<size_t, size_t> Ac_automata::match_pos_in_pattern(const std::string& pattern) {
    // return <start_pos_in_pat, length>

    std::vector<size_t> matched_min_factors_pos_in_text;
    Node* curr_node = root.get();
    for (size_t i = 0; i < pattern.size(); i++) {
        char c = pattern[i];
        
        // find the next node
        while (curr_node != nullptr && curr_node->children.find(c) == curr_node->children.end()) {
            curr_node = curr_node->fail;
        }

        if (curr_node == nullptr) {
            curr_node = root.get();
            continue;
        }

        curr_node = curr_node->children[c].get();
        
        if (curr_node->is_min_factor) {
            return std::make_pair(i - curr_node->height + 1, curr_node->height);
        }
    }
    return std::make_pair(-1, -1);
}

// output in binary format
void Ac_automata::serialize(std::string& filename) {
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
        for (const auto& p : curr_node->children) {
            s.push(p.second.get());
        }
    }

    // serialize
    index_type node_num = node_id.size();
    out.write(reinterpret_cast<char*>(&node_num), sizeof(node_num));
    for (const auto& p: node_id) {
        Node* node = p.first;

        index_type id = p.second;
        bool is_min_factor = node->is_min_factor; // optimize, done
        uint16_t height = node->height; // bound by lambda, 2^11
        index_type fail_link_id;
        uint8_t children_num = node->children.size(); // optimize, done
        uint8_t positions_num = node->positions.size(); // bound by tau, 2^5

        // combine is_min_factor and children_num into one byte
        uint8_t is_min_factor_children_num = (is_min_factor << 7) | children_num;
        // combine height and positions_num into one short
        uint16_t height_positions_num = (height << 5) | positions_num;
        
        // the fail link
        if (node->fail == nullptr) {
            fail_link_id = -1;
        } else {
            fail_link_id = node_id[node->fail];
        }

        // write primitive vars
        out.write(reinterpret_cast<char*>(&id), sizeof(id));
        // out.write(reinterpret_cast<char*>(&is_min_factor), sizeof(is_min_factor));
        // out.write(reinterpret_cast<char*>(&height), sizeof(height));
        out.write(reinterpret_cast<char*>(&fail_link_id), sizeof(fail_link_id));
        // out.write(reinterpret_cast<char*>(&children_num), sizeof(children_num));
        // out.write(reinterpret_cast<char*>(&positions_num), sizeof(positions_num)); // TODO optimize, fewer than tau

        // write combined vars
        out.write(reinterpret_cast<char*>(&is_min_factor_children_num), sizeof(is_min_factor_children_num));
        out.write(reinterpret_cast<char*>(&height_positions_num), sizeof(height_positions_num));

        for (const auto& p : node->children) {
            char c = p.first;
            index_type child_id = node_id[p.second.get()];
            out.write(reinterpret_cast<char*>(&c), sizeof(c)); // TODO optimize
            out.write(reinterpret_cast<char*>(&child_id), sizeof(child_id));
        }

        for (index_type pos : node->positions) {
            out.write(reinterpret_cast<char*>(&pos), sizeof(pos));
        }
    }
}

// read in binary format
void Ac_automata::load(std::string& filename) {
    using index_type = uint32_t;
    std::ifstream in(filename, std::ios::binary);

    // save id_node mapping, id_fail_link mapping, and id_children mapping, for first pass
    std::unordered_map<index_type, Node*> id_node;
    std::unordered_map<index_type, index_type> id_fail_link;
    std::unordered_map<index_type, std::vector<std::pair<char, index_type>>> id_children;

    // deserialize
    index_type node_num;
    in.read(reinterpret_cast<char*>(&node_num), sizeof(node_num));

    // first pass
    for (index_type i = 0; i < node_num; i++) {
        index_type id;
        bool is_min_factor; // improve
        uint16_t height; // bound by lambda
        index_type fail_link_id;
        uint8_t children_num; // fewer than sigma
        uint8_t positions_num; // bound by tau

        // combined vars
        uint8_t is_min_factor_children_num;
        uint16_t height_positions_num;

        // read primitive vars
        in.read(reinterpret_cast<char*>(&id), sizeof(id));
        // in.read(reinterpret_cast<char*>(&is_min_factor), sizeof(is_min_factor));
        // in.read(reinterpret_cast<char*>(&height), sizeof(height));
        in.read(reinterpret_cast<char*>(&fail_link_id), sizeof(fail_link_id));
        // in.read(reinterpret_cast<char*>(&children_num), sizeof(children_num));
        // in.read(reinterpret_cast<char*>(&positions_num), sizeof(positions_num));

        // read combined vars
        in.read(reinterpret_cast<char*>(&is_min_factor_children_num), sizeof(is_min_factor_children_num));
        is_min_factor = is_min_factor_children_num >> 7;
        children_num = is_min_factor_children_num & 0x7f;
        in.read(reinterpret_cast<char*>(&height_positions_num), sizeof(height_positions_num));
        height = height_positions_num >> 5;
        positions_num = height_positions_num & 0x001f;

        id_node[id] = new Node();
        id_node[id]->is_min_factor = is_min_factor;
        id_node[id]->height = height;
        id_fail_link[id] = fail_link_id;
        for (index_type j = 0; j < children_num; j++) {
            char c;
            index_type child_id;
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            in.read(reinterpret_cast<char*>(&child_id), sizeof(child_id));
            id_children[id].emplace_back(c, child_id);
        }
        for (index_type j = 0; j < positions_num; j++) {
            index_type pos;
            in.read(reinterpret_cast<char*>(&pos), sizeof(pos));
            id_node[id]->positions.emplace_back(pos);
        }
    }

    // second pass
    for (const auto& [id, node] : id_node) {
        if (id_fail_link[id] == -1) {
            node->fail = nullptr;
        } else {
            node->fail = id_node[id_fail_link[id]];
        }
        for (const auto& [c, child_id] : id_children[id]) {
            node->children[c] = std::unique_ptr<Node>(id_node[child_id]);
        }
    }

    root = std::unique_ptr<Node>(id_node[0]);
}