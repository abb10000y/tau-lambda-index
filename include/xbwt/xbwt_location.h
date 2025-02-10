#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <stack>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/int_vector.hpp>
#include <cstddef>
#include <memory>
#include "../symbol_table/symbol_bucket_offsets.h"
#include "../symbol_table/symbol_table.h"
#include <sdsl/qsufsort.hpp>
#include <chrono>

class XBWT_location
{
public:
    typedef std::pair<size_t, size_t> xbwt_range; // [l, r)

private:
    sdsl::wt_rlmn
    <
        sdsl::sd_vector<>,
        typename sdsl::sd_vector<>::rank_1_type,
        typename sdsl::sd_vector<>::select_1_type,
        sdsl::wt_ap<>
    > L;
    //sdsl::wt_ap<> L_wtap; // TODO: to be deleted
    sdsl::bit_vector P_bv;  // 1 for '(', 0 for ')'
    sdsl::bp_support_sada<> P;
    sdsl::bit_vector last;
    sdsl::bit_vector::rank_1_type last_rank;
    sdsl::bit_vector::select_1_type last_select;
    SymbolBucketOffsets C;
    uint64_t alphabet_size;
    uint64_t t_offset = static_cast<uint64_t>(1); // terminate symbols
    // sdsl::sd_vector<> singleCharRunCnt_bits;
    // sdsl::sd_vector<>::select_1_type singleCharRunCnt;
    // sdsl::sd_vector<> subFactorCnt_bits;
    // sdsl::sd_vector<>::select_1_type subFactorCnt;

    sdsl::int_vector<> locations, locations_offset, locations_length, next_leaf_rank; // leaves (rank only consider t_symbol)
    sdsl::int_vector<> leftest_child_rank, rightest_child_rank; // internal nodes (rank exclude t_symbol)
    SymbolTable symbol_table_;

    bool DownwardNavigation(xbwt_range& L_range, size_t c) const; // return found or not
    bool UpwardNavigation(size_t& l, size_t& r) const; // return found or not
    std::pair<size_t, size_t> pi_range_to_L_range(xbwt_range& pi_range) const; // pi range [l, r) to L/last range [l, r)
    std::pair<size_t, size_t> L_range_to_pi_range(xbwt_range& L_range, size_t c) const; // L/last range [l, r) to pi range [l, r) (starting with 'c')
    inline bool IsNotEmptyRange (xbwt_range const &range) const {
        return (std::get<0>(range) < std::get<1>(range));
    }

    void dfs(xbwt_range L_range, size_t &prev_leaf_rank, size_t prev_internal_node_rank);

    /*
    void print(xbwt_range range, bool isL) { // TODO: to be deleted
        if (isL) { std::cout << "[L]\t" ;}
        else { std::cout << "[pi]\t"; }
        std::cout << "l, r: " << std::get<0>(range) << ", " << std::get<1>(range) << "\n";
    }
    void processTime(std::chrono::steady_clock::time_point &t1, std::chrono::steady_clock::time_point &t2, std::string s); // TODO: to be deleted or integerated to "utility.h"
    */
public:
    XBWT_location();
    ~XBWT_location(){};

    void insert(sdsl::int_vector<> text);
    uint64_t match(sdsl::int_vector<>::iterator text_begin, sdsl::int_vector<>::iterator text_end);
    uint64_t nodeHeight(size_t l, size_t r);
    std::pair<size_t, size_t> match_pos_in_pattern(const sdsl::int_vector<>& pattern);
    bool match_if_exist(const sdsl::int_vector<>& pattern);
    void failureLink(size_t& l, size_t& r);
    void Serialize (std::ostream &out);
    //void Serialize_Partition (std::string const &path); // TODO: to be deleted
    void Load (std::istream &in);
    //void Load_Partition (std::string const &path); // TODO: to be deleted

    size_t getP_size() {return P.size();}
    std::vector<size_t> getC() {
        std::vector<size_t> ls;
        for (size_t i = 0; i < alphabet_size; i++)
            ls.push_back(C[i]);
        return ls;
    }
    std::vector<bool> getLast() {
        std::vector<bool> ls;
        for (auto i : last) ls.push_back(i);
        return ls;
    }
    std::vector<size_t> getL() {
        std::vector<size_t> ls;
        for (auto i : L) ls.push_back(i);
        return ls;
    }
    uint64_t getGrammarNumber() { return L.rank(L.size(), 1); }
    uint16_t getAlphabetCnt(uint64_t c) { return C[c+1] - C[c]; }
    void getPiPath(size_t i, sdsl::int_vector<>& path) const; // node-to-root path for the i-th '$'
    // template <typename Iterator, typename Range>
    // void suffixRangeQuery(Iterator rbegin, Iterator rend, Range &grammar_range, SymbolTable& runLength_symbol_table, uint64_t maxRunLength);
    // template <typename Iterator, typename Range>
    // void prefixRangeQuery(Iterator begin, Iterator end, Range &grammar_range, SymbolTable& runLength_symbol_table, uint64_t maxRunLength);
    // template <typename Iterator>
    // uint64_t infixMatch(Iterator begin, Iterator end); // assume the input is in lex-order
    // void gen_subFactorCnt(sdsl::int_vector<> &lex_sl_rl_int_text, SymbolTable &symbol_table, uint64_t maxRunLength);
    // template <typename Iterator>
    // uint64_t subFactorSuffixRangeQuery(Iterator rbegin, Iterator rend, SymbolTable& runLength_symbol_table, uint64_t maxRunLength);
    // void singleCharacterRunQuery(uint64_t rl_int_symbol, SymbolTable& runLength_symbol_table, std::pair<uint64_t, uint64_t> &result, uint64_t maxRunLength); // {cnt, len}
    friend std::ostream& operator<< (std::ostream &out, XBWT_location const &xbwt);
    void insert_locations(std::vector<std::vector<size_t>> &locations_tmp, size_t locations_cnt);
};

void XBWT_location::dfs(xbwt_range L_range, size_t &prev_leaf_rank, size_t prev_internal_node_rank) {
    size_t leaf_rank = L.rank(std::get<1>(L_range), 1);
    for (size_t i = std::get<0>(L_range); i < std::get<1>(L_range); i++) {
        size_t c = L[i];
        if (c == 1) {
            next_leaf_rank[prev_leaf_rank] = leaf_rank;
            prev_leaf_rank = leaf_rank;
            if (std::get<0>(L_range) > 0) { // current node is not the root
                leftest_child_rank[prev_internal_node_rank] = leaf_rank;
                rightest_child_rank[prev_internal_node_rank] = leaf_rank;
            }
        } else {
            xbwt_range cur_range = L_range;
            DownwardNavigation(cur_range, c);
            dfs(cur_range, prev_leaf_rank, i - leaf_rank);
            if (std::get<0>(L_range) > 0) { // current node is node the root
                if (i == std::get<0>(L_range)) { leftest_child_rank[prev_internal_node_rank] = leftest_child_rank[i - leaf_rank]; }
                rightest_child_rank[prev_internal_node_rank] = rightest_child_rank[i - leaf_rank];
            }
        }
    }
}

void XBWT_location::insert_locations(std::vector<std::vector<size_t>> &locations_tmp, size_t locations_cnt) {
    locations_offset.resize(locations_tmp.size());
    locations_length.resize(locations_tmp.size());
    next_leaf_rank.resize(locations_tmp.size());
    locations.resize(locations_cnt);
    for (size_t i = 0, j = 0, k = 0, cnt = 0, end = locations_tmp.size(); i < end; i++, j++) {
        locations_offset[j] = cnt;
        locations_length[j] = locations_tmp[i].size();
        cnt += locations_tmp[i].size();
        for (auto loc : locations_tmp[i]) { locations[k++] = loc; }
    }

    leftest_child_rank.resize(last.size() - locations_tmp.size());
    rightest_child_rank.resize(last.size() - locations_tmp.size());
    xbwt_range L_range {0, last_select(1) + 1};
    size_t prev_leaf_rank = 0;
    dfs(L_range, prev_leaf_rank, 0);

    sdsl::util::bit_compress(locations_offset);
    sdsl::util::bit_compress(locations_length);
    sdsl::util::bit_compress(next_leaf_rank);
    sdsl::util::bit_compress(leftest_child_rank);
    sdsl::util::bit_compress(rightest_child_rank);
}

XBWT_location::XBWT_location() {}

void XBWT_location::failureLink(size_t& l, size_t& r) {
    /*
    if (i == 0) return 0; // empty string
    size_t k = P.select(i+1);
    size_t j = P.enclose(k);
    size_t l = P.rank(j)-1;
    */
    if (l == 0) return; // empty string
    size_t k = P.select(last_rank(l) + 1);
    size_t j = P.enclose(k);
    if (P.rank(j) - 1 == 0) { l = 0; }
    else { l = last_select(P.rank(j) - 1) + 1; }
    r = last_select(last_rank(l) + 1) + 1;
}

/*
void XBWT_location::processTime(std::chrono::steady_clock::time_point &t1, std::chrono::steady_clock::time_point &t2, std::string s) {
    t2 = std::chrono::steady_clock::now();
        // std::cout << s << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << " ns"; // "comsum time (ns):;
        std::cout << s << std::chrono::duration_cast<std::chrono::minutes>(t2 - t1).count() << " mins"; // "comsum time (mins):;
        std::cout << "\n";
    t1 = std::chrono::steady_clock::now();
}
*/

void XBWT_location::insert(sdsl::int_vector<> text) {
    // appending '0' to text is necessary for the s.sort(), however we don't want '0' appearing in xbwt
    // hence we will skip 'the 1st element in sa' and 'the last element in text'
    if (text[text.size()-1] != 0) { sdsl::append_zero_symbol(text); }
    symbol_table_ = decltype(symbol_table_)(text);
    for (size_t i = 0; i < text.size(); i++) { text[i] = symbol_table_[text[i]]; }
    alphabet_size = symbol_table_.GetAlphabetSize();
    sdsl::int_vector<> sa;
    sdsl::qsufsort::construct_sa(sa, text);

    // MR, modified Kasai alg.
    std::vector<size_t> RCP;
    std::vector<bool> MR;
    size_t len = text.size() - 1; // negelect the '0'
    sdsl::int_vector<> inv_sa(len);
    MR.assign(len, false);
    RCP.assign(len, 0);
    MR[0] = 1;
    for (size_t i = 0; i < len; i++) inv_sa[sa[i+1]] = i; // +1 for skipping the '0'
    for (size_t i = 0, lcp = 0, rcp = 0; i < len-1; ++i) {
        int j = sa[inv_sa[i] - 1 + 1]; // -1 for the previous position; +1 for skipping the '0'
        while (i + lcp < len && j + lcp < len && text[i + lcp] == text[j + lcp])
            ++lcp;
        while (i + rcp < len && j + rcp < len && text[i + rcp] != t_offset && text[i + rcp] == text[j + rcp])
            ++rcp;
        MR[inv_sa[i]] = lcp <= rcp;
        RCP[inv_sa[i]] = rcp;
        if (lcp > 0) --lcp;
        if (rcp > 0) --rcp;
    }

    // Construct L, Last, and P from MR

    /*
    // TODO: a better alg. for constructing P
    auto isPrefix = [&](int a, int b) {
        int i = a, j = b;
        while (text[i] != t_offset && text[j] != t_offset && text[i] == text[j]) {
            i++;
            j++;
        }
        return text[i] == t_offset;
    };
    */
//processTime(t1, t2, "chk7_2: ");

    std::vector<uint64_t> L_vector;
    std::vector<bool> last_vector, P_vector;
    sdsl::int_vector<> C_vector;
    // std::stack<size_t> stk; // store the first index of pi
    std::stack<std::pair<size_t, size_t>> stk;
    MR.push_back(1);
    for (size_t i = 0, prev = 0; i <= len; i++) {
        if (MR[i] == 1) {
            // L and Last
            std::set<uint64_t> st;
            for (size_t j = prev; j < i; j++) {
                size_t idx = sa[j + 1];  // +1 for skipping the '0'
                st.insert(text[idx > 0 ? idx-1 : len-1]); // BWT
            }
            for (auto itr = st.begin(), e = --st.end(); itr != st.end(); itr++) {
                L_vector.push_back(*itr);
                last_vector.push_back(itr == e);
            }
            prev = i;

            // P
            if (i < len) {
                /*
                while (!stk.empty() && !isPrefix(stk.top(), sa[i + 1])) {  // +1 for skipping the '0'
                    P_vector.push_back(0);
                    stk.pop();
                }
                P_vector.push_back(1);
                stk.push(sa[i + 1]);  // +1 for skipping the '0'
                */
                size_t l = RCP[i];
                // size > 1 for the first element is empty string (less than 1)
                // which should be the prefix of all substrings (always in stack)
                while (stk.size() > 1 && stk.top().second >= l) { 
                    P_vector.push_back(0);
                    stk.pop();
                }
                if (stk.size() > 1 && text[sa[stk.top().first + 1] + l] != t_offset) {
                    l = stk.top().second;
                    P_vector.push_back(0);
                    stk.pop();
                }
                P_vector.push_back(1);
                stk.push({i, l});
            }
        }
    }
    while (!stk.empty()){
        P_vector.push_back(0);
        stk.pop();
    }
    MR.pop_back();
//processTime(t1, t2, "chk7_3: ");

    // to succint data structure
    last.resize(last_vector.size());
    for (size_t i = 0; i < last_vector.size(); i++) last[i] = last_vector[i];
    last_rank = decltype(last_rank)(&last);
    last_select = decltype(last_select)(&last);
//processTime(t1, t2, "chk7_4: ");

    sdsl::int_vector<> L_buffer;
    L_buffer.width(symbol_table_.GetEffectiveAlphabetWidth());
    L_buffer.resize(L_vector.size());
    for (size_t i = 0; i < L_vector.size(); i++) L_buffer[i] = L_vector[i];
    sdsl::construct_im(L, L_buffer);
    //sdsl::construct_im(L_wtap, L_buffer); // TODO: to be deleted
//processTime(t1, t2, "chk7_5: ");

    sdsl::bit_vector bv_tmp (P_vector.size(), 0);
    for (size_t i = 0; i < P_vector.size(); i++) bv_tmp[i] = P_vector[i];
    P_bv = bv_tmp;
    sdsl::bp_support_sada<> tmp(&P_bv);
    P = tmp;
//processTime(t1, t2, "chk7_6: ");
    C_vector.resize(alphabet_size + 1); // +1 for counting the # of last alphabet
    for (size_t i = 1, nodeCnt = 1; i <= alphabet_size; i++) {
        size_t cnt = L.rank(L.size(), i-1);
        nodeCnt += cnt;
        C_vector[i] = nodeCnt;
    }
    C = decltype(C)(C_vector);
//processTime(t1, t2, "chk7_7: ");
}

// void XBWT_location::gen_subFactorCnt(sdsl::int_vector<> &lex_sl_rl_int_text, SymbolTable &symbol_table, uint64_t maxRunLength) {
//     sdsl::int_vector<> subFactorCnt_bits_tmp(last_rank(last.size()), 0);
//     sdsl::int_vector<> singleCharRunCnt_bits_tmp(last_rank(last.size()), 0);
//     auto rbegin {std::prev(lex_sl_rl_int_text.end(), 3)}; // skipping the last '1''0'
//     auto rend {std::prev(lex_sl_rl_int_text.begin())};
//     auto rfirst {rbegin};
//     auto rlast {rfirst};
//     while (rlast != rend) {
//         while (rlast != rend && *rlast != t_offset) { rlast = std::prev(rlast); }
//         auto it {rfirst};
//         // match function from right to left
//         xbwt_range L_range = {0, last_select(1) + 1};
//         while (it != rlast) {
//             size_t c = *it--;
//             if (!DownwardNavigation(L_range, c)) {
//                 throw std::runtime_error("matching in gen_subFactorCnt() should never fail");
//             }
//             subFactorCnt_bits_tmp[last_rank(std::get<1>(L_range)) - 1] += 1; // -1 for 0-index
//             singleCharRunCnt_bits_tmp[last_rank(std::get<1>(L_range)) - 1] += 1; // -1 for 0-index
//         }
//         if (rlast != rend) {
//             rlast = std::prev(rlast);
//             rfirst = rlast;
//         }
//     }
//     //std::cout << "symbol_table: " << symbol_table << "\n";
//     for (size_t i = 2; i < alphabet_size; i++) {
//         size_t mul = symbol_table.get_run_length(i, maxRunLength);
//         size_t begin = C[i] - C[2] + 1, end = C[i+1] - C[2] + 1;
//         for (size_t j = begin; j < end; j++) {
//             singleCharRunCnt_bits_tmp[j] = mul * singleCharRunCnt_bits_tmp[j];
//         }
//     }
//     for (size_t i = 1; i < subFactorCnt_bits_tmp.size(); i++) {
//         subFactorCnt_bits_tmp[i] += subFactorCnt_bits_tmp[i-1];
//         singleCharRunCnt_bits_tmp[i] += singleCharRunCnt_bits_tmp[i-1];
//     }
//     //std::cout << "lex_sl_rl_int_text: " << lex_sl_rl_int_text << "\n";
//     //std::cout << "subFactorCnt_bits_tmp: " << subFactorCnt_bits_tmp << "\n";
//     //std::cout << "singleCharRunCnt_bits_tmp: " << singleCharRunCnt_bits_tmp << "\n";
//     subFactorCnt_bits = decltype(subFactorCnt_bits)(std::begin(subFactorCnt_bits_tmp), std::end(subFactorCnt_bits_tmp));
//     subFactorCnt.set_vector(&subFactorCnt_bits);
//     singleCharRunCnt_bits = decltype(singleCharRunCnt_bits)(std::begin(singleCharRunCnt_bits_tmp), std::end(singleCharRunCnt_bits_tmp));
//     singleCharRunCnt.set_vector(&singleCharRunCnt_bits);
// }

void XBWT_location::Serialize(std::ostream &out) {
    sdsl::write_member(alphabet_size, out);
    L.serialize(out);
    P_bv.serialize(out);
    P.serialize(out);
    last.serialize(out);
    last_rank.serialize(out);
    last_select.serialize(out);
    C.Serialize(out);
    // subFactorCnt_bits.serialize(out);
    // subFactorCnt.serialize(out);
    // singleCharRunCnt_bits.serialize(out);
    // singleCharRunCnt.serialize(out);
    locations_offset.serialize(out);
    locations_length.serialize(out);
    next_leaf_rank.serialize(out);
    leftest_child_rank.serialize(out);
    rightest_child_rank.serialize(out);
    symbol_table_.Serialize(out);
}

void XBWT_location::Load(std::istream &in) {
    sdsl::read_member(alphabet_size, in);
    L.load(in);
    P_bv.load(in);
    P.load(in, &P_bv);
    last.load(in);
    last_rank.load(in, &last);
    last_select.load(in, &last);
    C.Load(in);
    // subFactorCnt_bits.load(in);
    // subFactorCnt.load(in, &subFactorCnt_bits);
    // singleCharRunCnt_bits.load(in);
    // singleCharRunCnt.load(in, &singleCharRunCnt_bits);
    locations_offset.load(in);
    locations_length.load(in);
    next_leaf_rank.load(in);
    leftest_child_rank.load(in);
    rightest_child_rank.load(in);
    symbol_table_.Load(in);
}

// match [begin, ..., end)
uint64_t XBWT_location::match(sdsl::int_vector<>::iterator text_begin, sdsl::int_vector<>::iterator text_end) {
    if (last_rank(last_rank.size()) == 0) return 0; // empty
    xbwt_range L_range = {0, last_select(1) + 1};
    while (text_begin != text_end) {
        size_t c = symbol_table_[*text_begin++]; ;
        if (!DownwardNavigation(L_range, c)) { return 0; }
    }
    size_t grammar_below_l = L.rank(std::get<0>(L_range), t_offset);
    size_t grammar_below_r = L.rank(std::get<1>(L_range), t_offset);
    if (grammar_below_l == grammar_below_r) { return 0; } // no grammar within [l, r)
    else { return grammar_below_r; }
}

bool XBWT_location::DownwardNavigation(xbwt_range& L_range, size_t c) const{
    xbwt_range pi_range = L_range_to_pi_range(L_range, c);
    if (!IsNotEmptyRange(pi_range)) { return false; }
    L_range = pi_range_to_L_range(pi_range);
    return true;
}

bool XBWT_location::UpwardNavigation(size_t& l, size_t& r) const{
    if (l == 0) return false; // root (empty string)
    size_t nodeCnt = last_rank(r) - 1; // -1 for 0-index 
    // size_t c = C.getSymbol(nodeCnt + 1 + C[2] - 1); // +1 for rank (rightmost excluded)
    size_t c = C.getSymbol(nodeCnt+ C[2]);
    //size_t c_cnt = nodeCnt - (C[c] - C[2] + 1) + 1; // 0-index
    size_t c_cnt = nodeCnt - C[c] + C[2];
    if (last_rank(L.select(c_cnt, c)) == 0) { l = 0; }
    else { l = last_select(last_rank(L.select(c_cnt, c))) + 1; }
    // r = last_select(last_rank(L.select(cnt, c)) + 1) + 1;
    r = last_select(last_rank(l) + 1) + 1;
    return true;
}

uint64_t XBWT_location::nodeHeight(size_t l, size_t r) {
    uint64_t height = 0;
    while (UpwardNavigation(l, r)) { height++; }
    return height;
}

void XBWT_location::getPiPath(size_t i, sdsl::int_vector<>& path) const {
    if (i == 0 || i > L.rank(L.size(), t_offset)) { return; }
    std::vector<uint64_t> tmp;
    size_t l = last_select(last_rank(L.select(i, t_offset))) + 1;
    size_t r = last_select(last_rank(l) + 1) + 1;
    while (l > 0) {
        tmp.push_back(C.getSymbol(last_rank(r) - 1 + 1 + C[2] - 1));
        UpwardNavigation(l, r);
    }
    path.resize(tmp.size());
    for (size_t i = 0; i < tmp.size(); i++) {
        path[i] = tmp[i];
    }
}

// void XBWT_location::singleCharacterRunQuery(uint64_t rl_int_symbol, SymbolTable& runLength_symbol_table, std::pair<uint64_t, uint64_t> &result, uint64_t maxRunLength) {
//     uint64_t l_rl_int_symbol = rl_int_symbol, r_rl_int_symbol = runLength_symbol_table.max_of_this_alphabet(l_rl_int_symbol, maxRunLength); // [l, r)
//     xbwt_range pi_range = {
//         C[l_rl_int_symbol] - C[2] + 1,
//         C[r_rl_int_symbol] - C[2] + 1
//     };
//     // std::cout << "l/r_rl_int_symbol: " << l_rl_int_symbol << ", " << r_rl_int_symbol << "\n";
//     // std::cout << "pi_range: " << std::get<0>(pi_range) << ", " << std::get<1>(pi_range) << "\n";
//     // {cnt * len, cnt}
//     result = {
//         singleCharRunCnt(std::get<1>(pi_range)) - singleCharRunCnt(std::get<0>(pi_range)),
//         subFactorCnt(std::get<1>(pi_range)) - subFactorCnt(std::get<0>(pi_range))
//     };
//     // std::cout << subFactorCnt(std::get<0>(pi_range)) << ", ";
//     // std::cout << subFactorCnt(std::get<1>(pi_range)) << ", ";
//     // std::cout << singleCharRunCnt(std::get<0>(pi_range)) << ", ";
//     // std::cout << singleCharRunCnt(std::get<1>(pi_range)) << "\n";
//     //return singleCharRunCnt(C[r_rl_int_symbol] - C[2] + 1) - singleCharRunCnt(C[l_rl_int_symbol] - C[2] + 1);
// }

// template <typename Iterator>
// uint64_t XBWT_location::subFactorSuffixRangeQuery(Iterator rbegin, Iterator rend, SymbolTable& runLength_symbol_table, uint64_t maxRunLength) {
//     uint64_t l_rl_int_symbol = *rbegin, r_rl_int_symbol = runLength_symbol_table.max_of_this_alphabet(l_rl_int_symbol, maxRunLength); // [l, r)
//     xbwt_range L_range = {
//         last_select(C[l_rl_int_symbol] - C[2] + 1) + 1,
//         last_select(C[r_rl_int_symbol] - C[2] + 1) + 1
//     };
//     // std::cout << "L_range: " << std::get<0>(L_range) << ", " << std::get<1>(L_range) << "\n";
//     size_t cnt = 0;
//     rend = std::next(rend); // range-query until the last character run
//     auto it {std::prev(rbegin)};
//     while (it != rend) {
//         //print(L_range, true);
//         if (!DownwardNavigation(L_range, *it)) { return cnt; }
//         it = std::prev(it);
//     }

//     l_rl_int_symbol = *rend, r_rl_int_symbol = runLength_symbol_table.max_of_this_alphabet(l_rl_int_symbol, maxRunLength); // [l, r)
//     //print(L_range, true);
//     // the last character run
//     for (size_t rl_int_symbol = l_rl_int_symbol; rl_int_symbol < r_rl_int_symbol; rl_int_symbol++) {
//         xbwt_range tmp_L_range {L_range};
//         //std::cout << "len: " << len << "\n";
//         xbwt_range pi_range = L_range_to_pi_range(tmp_L_range, rl_int_symbol);
//         //print(pi_range, false);
//         if (!IsNotEmptyRange(pi_range)) { continue; }
//         cnt += subFactorCnt(std::get<1>(pi_range)) - subFactorCnt(std::get<0>(pi_range));
//     }
//     return cnt;
// }

// template <typename Iterator, typename Range>
// void XBWT_location::suffixRangeQuery(Iterator rbegin, Iterator rend, Range &grammar_range, SymbolTable& runLength_symbol_table, uint64_t maxRunLength) {
//     uint64_t l_rl_int_symbol = *rbegin, r_rl_int_symbol = runLength_symbol_table.max_of_this_alphabet(l_rl_int_symbol, maxRunLength); // [l, r)
//     xbwt_range L_range = {
//         last_select(C[l_rl_int_symbol] - C[2] + 1) + 1,
//         last_select(C[r_rl_int_symbol] - C[2] + 1) + 1
//     };
//     // std::cout << "L_range: " << std::get<0>(L_range) << ", " << std::get<1>(L_range) << "\n";
//     // xbwt_range L_range = {last_select(C[*rbegin] - C[2] + 1) + 1, last_select(C[(*rbegin) + 1] - C[2] + 1)};
//     auto it {std::prev(rbegin)};
//     while (it != rend) {
//         if (!DownwardNavigation(L_range, *it)) { 
//             grammar_range = {1, 0}; 
//             return;
//         }
//         it = std::prev(it);
//     }
//     // std::cout << "L_range: " << std::get<0>(L_range) << ", " << std::get<1>(L_range) << "\n";
//     // std::cout << "l, r: " << L.rank(std::get<0>(L_range), t_offset) + 1 << ", " << L.rank(std::get<1>(L_range), t_offset) << "\n";
//     grammar_range = {L.rank(std::get<0>(L_range), t_offset) + 1, L.rank(std::get<1>(L_range), t_offset)};
//     if (std::get<0>(grammar_range) > std::get<1>(grammar_range)) {
//         grammar_range = {1, 0};
//     }
//     return;
// }

// template <typename Iterator, typename Range>
// void XBWT_location::prefixRangeQuery(Iterator begin, Iterator end, Range &grammar_range, SymbolTable& runLength_symbol_table, uint64_t maxRunLength) {
//     uint64_t l_rl_int_symbol = *begin, r_rl_int_symbol = runLength_symbol_table.max_of_this_alphabet(l_rl_int_symbol, maxRunLength); // [l, r)
//     xbwt_range L_range = {
//         last_select(C[l_rl_int_symbol] - C[2] + 1) + 1,
//         last_select(C[r_rl_int_symbol] - C[2] + 1) + 1
//     };
//     //xbwt_range L_range = {last_select(C[*begin] - C[2] + 1) + 1, last_select(C[(*begin) + 1] - C[2] + 1)};
//     auto it {std::next(begin)};
//     while (it != end) {
//         if (!DownwardNavigation(L_range, *it)) { 
//             grammar_range = {1, 0}; 
//             return;
//         }
//         it = std::next(it);
//     }
//     grammar_range = {L.rank(std::get<0>(L_range), t_offset) + 1, L.rank(std::get<1>(L_range), t_offset)};
//     if (std::get<0>(grammar_range) > std::get<1>(grammar_range)) {
//         grammar_range = {1, 0};
//     }
//     return;
// }

// template <typename Iterator>
// uint64_t XBWT_location::infixMatch(Iterator begin, Iterator end) {
//     xbwt_range L_range = {0, last_select(1) + 1};
//     auto it {begin};
//     while (it != end) {
//         if (!DownwardNavigation(L_range, *it)) { return 0; }
//         it = std::next(it);
//     }
//     size_t grammar_below_l = L.rank(std::get<0>(L_range), t_offset);
//     size_t grammar_below_r = L.rank(std::get<1>(L_range), t_offset);
//     if (grammar_below_l == grammar_below_r) { return 0; } // no grammar within [l, r)
//     else { return grammar_below_r; }
// }

std::ostream& operator<< (std::ostream &out, XBWT_location const &xbwt) {
    size_t n = xbwt.last.size();
    out << "idx\tlast\tL\tpi\n";
    for (size_t i = 0; i < n; i++) {
        out << i << "\t";
        out << xbwt.last[i] << "\t";
        out << xbwt.L[i] << "\t";
        if (i > 0 && xbwt.last[i-1] == 1) {
            size_t l = xbwt.last_select(xbwt.last_rank(i)) + 1;
            size_t r = xbwt.last_select(xbwt.last_rank(l) + 1) + 1;
            while (l > 0) {
                out << xbwt.C.getSymbol(xbwt.last_rank(r) - 1 + 1 + xbwt.C[2] - 1) << " ";
                xbwt.UpwardNavigation(l, r);
            }
        }
        out << std::endl;
    }
    return out;
}

std::pair<size_t, size_t> XBWT_location::pi_range_to_L_range(xbwt_range& pi_range) const {
    // size_t L_left = last_select(std::get<0>(pi_range)) + 1, L_right = last_select(std::get<1>(pi_range)) + 1;
    return {
        last_select(std::get<0>(pi_range)) + 1,
        last_select(std::get<1>(pi_range)) + 1
    };
}

std::pair<size_t, size_t> XBWT_location::L_range_to_pi_range(xbwt_range& L_range, size_t c) const {
    // -(C[2]-1) shifting for skipping the C[1]('$'); (L.rank(l,c)+1) to count the i-th c in pi; -1 for the 0-index in pi
    // size_t pi_left = C[c] - (C[2] - 1) + (L.rank(l, c) + 1) - 1;
    // size_t pi_left = C[c] - C[2] + L.rank(std::get<0>(L_range), c) + 1, pi_right = C[c] - C[2] + L.rank(std::get<1>(L_range), c) + 1;
    return {
        C[c] - C[2] + L.rank(std::get<0>(L_range), c) + 1,
        C[c] - C[2] + L.rank(std::get<1>(L_range), c) + 1
    };
}

/*
void XBWT_location::Serialize_Partition (std::string const &path) {
    {
        std::fstream out {path + "_L_wt_rlmn", std::ios_base::out | std::ios_base::trunc};
        sdsl::write_member(alphabet_size, out);
        L.serialize(out);
        out.close();
    }
    {
        std::fstream out {path + "_L_wt_ap", std::ios_base::out | std::ios_base::trunc};
        L_wt_ap.serialize(out);
        out.close();
    }
    {
        std::fstream out {path + "_P", std::ios_base::out | std::ios_base::trunc};
        P_bv.serialize(out);
        P.serialize(out);
        out.close();
    }
    {
        std::fstream out {path + "_last", std::ios_base::out | std::ios_base::trunc};
        last.serialize(out);
        last_rank.serialize(out);
        last_select.serialize(out);
        out.close();
    }
    {
        std::fstream out {path + "_C", std::ios_base::out | std::ios_base::trunc};
        C.Serialize(out);
        out.close();
    }
    {
        std::fstream out {path + "_subFactorCnt_bits", std::ios_base::out | std::ios_base::trunc};
        subFactorCnt_bits.serialize(out);
        subFactorCnt.serialize(out);
        out.close();
    }
}
*/

/*
void XBWT_location::Load_Partition (std::string const &path) {
    {
        std::fstream in {path + "_L_wt_rlmn", std::ios_base::in};
        sdsl::read_member(alphabet_size, in);
        L.load(in);
        in.close();
    }
    {
        std::fstream in {path + "_L_wt_ap", std::ios_base::in};
        L_wtap.load(in);
        in.close();
    }
    {
        std::fstream in {path + "_P", std::ios_base::in};
        P_bv.load(in);
        P.load(in, &P_bv);
        in.close();
    }
    {
        std::fstream in {path + "_last", std::ios_base::in};
        last.load(in);
        last_rank.load(in, &last);
        last_select.load(in, &last);
        in.close();
    }
    {
        std::fstream in {path + "_C", std::ios_base::in};
        C.Load(in);
        in.close();
    }
    {
        std::fstream in {path + "_subFactorCnt_bits", std::ios_base::in};
        subFactorCnt_bits.load(in);
        subFactorCnt.load(in, &subFactorCnt_bits);
        in.close();
    }
}
*/
std::pair<size_t, size_t> XBWT_location::match_pos_in_pattern(const sdsl::int_vector<>& pattern) {
    // uint64_t t_offset = static_cast<uint64_t>(1); // terminate symbols // TODO: change to global variable
    // sdsl::int_vector<> pattern_int;
    // pattern_int.resize(pattern.size());
    // for (size_t i = 0; i < pattern.size(); i++) { pattern_int[i] = pattern[i] + t_offset; }

    xbwt_range L_range = {0, last_select(1) + 1};
    //size_t l = 0, r = last_select(1) + 1; // rank will exclusive the right most position
    for (size_t i = 0; i < pattern.size(); i++) {
        size_t c = symbol_table_[pattern[i]];
        while (!DownwardNavigation(L_range, c) && std::get<0>(L_range) > 0) { failureLink(std::get<0>(L_range), std::get<1>(L_range)); }
        if (L[std::get<0>(L_range)] == 1) { 
            size_t len = nodeHeight(std::get<0>(L_range), std::get<1>(L_range));
            return std::make_pair(i - len + 1, len);
        }
    }
    return std::make_pair(-1, -1);
}

bool XBWT_location::match_if_exist(const sdsl::int_vector<>& pattern) {
    xbwt_range L_range = {0, last_select(1) + 1};
    bool result = false;
    //size_t l = 0, r = last_select(1) + 1; // rank will exclusive the right most position
    for (size_t i = 0; i < pattern.size(); i++) {
        size_t c = symbol_table_[pattern[i]];
        while (!DownwardNavigation(L_range, c) && std::get<0>(L_range) > 0) { failureLink(std::get<0>(L_range), std::get<1>(L_range)); }
        if (L[std::get<0>(L_range)] == 1) { 
            result = true;
            break;
        }
    }
    return result;
}