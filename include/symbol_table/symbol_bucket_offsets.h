#pragma once

#include "../../util/utility.h"
#include <utility>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

class SymbolBucketOffsets
{
public:

  SymbolBucketOffsets () = default;
  SymbolBucketOffsets (SymbolBucketOffsets const&) = delete;
  SymbolBucketOffsets (SymbolBucketOffsets&&);
  SymbolBucketOffsets (sdsl::int_vector<> const &offsets);
  SymbolBucketOffsets& operator= (SymbolBucketOffsets const&) = delete;
  SymbolBucketOffsets& operator= (SymbolBucketOffsets&&);

  void Swap (SymbolBucketOffsets&);

  uint64_t Serialize
  (
    std::ostream &out,
    std::shared_ptr<SpaceNode> parent = nullptr,
    std::string const name = ""
  );
  void Load (std::istream &in);

  uint64_t operator[] (uint64_t const symbol) const;
  
  uint64_t at (uint64_t const symbol) const;

  uint64_t getSymbol (uint64_t pos) const; // the largest symbol in [0, pos)

  friend std::ostream& operator<< (std::ostream &out, SymbolBucketOffsets const &symbol_bucket_offsets);

private:

  sdsl::sd_vector<> offset_bits_;
  sdsl::sd_vector<>::select_1_type offset_select_1_;
  sdsl::sd_vector<>::rank_1_type offset_rank_1_;

};