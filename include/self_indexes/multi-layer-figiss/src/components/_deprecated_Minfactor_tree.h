#include <sdsl/cst_sct3.hpp>
#include <string>

// utilize cst in sdsl to store the min_factors
class Minfactor_tree {
public:
    Minfactor_tree(std::string_view text, vector<pair<size_t, size_t>>& min_factors);
    ~Minfactor_tree();

    // given one pattern, return all the positions of occurences of minimal factors
    vector<pair<size_t, size_t>> search(std::string_view pattern);
private:
    constexpr static char TERMINATOR = '$'; // TODO better way to define this?

    sdsl::cst_sct3<> cst;
    std::string_view text;

    using cst_id_type = typename sdsl::cst_sct3<>::size_type;
    std::map<cst_id_type, std::vector<pair<size_t, size_t>>> id_to_pos; // Task: size_t to pair<size_t, size_t> !!!
    std::set<cst_id_type> begin_from_factor;
    std::set<cst_id_type> is_end;

    void build_cst(const std::map<std::string_view, std::vector<size_t>>& min_factors_map);
};

Minfactor_tree::Minfactor_tree(std::string_view text, vector<pair<size_t, size_t>>& min_factors) : text(text) {
    // remove duplicate minimal factors
    std::map<std::string_view, vector<size_t>> min_factors_map;

    for (auto p : min_factors) {
        std::string_view factor = text.substr(p.first, p.second - p.first + 1);
        if (min_factors_map.find(factor) == min_factors_map.end()) {
            min_factors_map[factor] = vector<size_t>();
        }
        min_factors_map[factor].push_back(p.first);
    }

    // build the CST
    build_cst(min_factors_map);
}

Minfactor_tree::~Minfactor_tree() {
}

void Minfactor_tree::build_cst(const std::map<std::string_view, std::vector<size_t>>& min_factors_map) {
    // concatenate all the minimal factors
    std::string min_factors_str = "";
    for (auto p : min_factors_map) {
        min_factors_str += p.first;
        min_factors_str += TERMINATOR;
    }

    // build the CST
    sdsl::construct_im(cst, min_factors_str, 1); // TODO what is 1?

    // set id_to_pos, begin_from_factor, is_end
    for (auto p : min_factors_map) {
        std::string_view factor = p.first;
        vector<size_t> positions = p.second;

        auto curr_node = cst.root();
        for (size_t i = 0; i < factor.size(); i++) {
            // traverse the CST
            curr_node = cst.select_child(curr_node, factor[i]);
            auto id = cst.id(curr_node);
            begin_from_factor.insert(id);

            // reach the end of the factor, record the positions
            if (i == factor.size() - 1) {
                is_end.insert(id);
                for (auto pos : positions) {
                    id_to_pos[id].push_back(pos);
                }
            }
        }
    }
}

// given one pattern, return all the positions of occurences of minimal factors
bool Minfactor_tree::search(std::string_view pattern, vector<pair<size_t, size_t>& occ_text, vector<pair<size_t, size_t>& occ_pattern) {
    auto curr_node = cst.root();
    for (int l, r = 0; r < pattern.size(); r++) {
        if (cst.child(curr_node, pattern[r]) == cst::END) {
            return false;
        }
    }
}
