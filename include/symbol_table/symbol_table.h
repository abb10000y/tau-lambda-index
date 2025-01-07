#pragma once

#include "../util/utility.h"
#include <utility>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>

class SymbolTable
{
public:

  SymbolTable () = default;
  SymbolTable (SymbolTable const&) = delete;
  SymbolTable (SymbolTable&&);
  SymbolTable (sdsl::int_vector<> const &byte_text);
  SymbolTable& operator= (SymbolTable const&) = delete;
  SymbolTable& operator= (SymbolTable &&);

  void Swap (SymbolTable&);
  
  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  inline uint8_t GetEffectiveAlphabetWidth () const
  {
    return effective_alphabet_width_;
  }

  // void SetAlphabetSize(uint64_t s) { alphabet_size = s; }
  uint64_t GetAlphabetSize() { return alphabet_size; }

  uint64_t operator[] (uint64_t const byte) const;

  uint64_t lower_bound(uint64_t int_symbol, uint64_t length, uint64_t maxLength) const;

  uint64_t max_of_this_alphabet(uint64_t rl_int_symbol, uint64_t maxLength) const;

  uint64_t symbol_to_byte(uint64_t const symbol) const { return alphabet_select_1_(symbol); }

  friend std::ostream& operator<< (std::ostream &out, SymbolTable const &symbol_table);

  uint64_t get_run_length(uint64_t c, uint64_t maxRunLength);

private:

  uint64_t effective_alphabet_width_;
  uint64_t alphabet_size;
  sdsl::bit_vector alphabet_bits_;
  sdsl::bit_vector::rank_1_type alphabet_rank_1_;
  sdsl::bit_vector::select_1_type alphabet_select_1_;

};