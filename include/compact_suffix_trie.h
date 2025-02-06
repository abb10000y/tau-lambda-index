#pragma once

#include "xbwt/xbwt.h"

class compact_suffix_trie {
private:
    XBWT *forward_xbwt {new XBWT()}, *reverse_xbwt {new XBWT()};
public:
    compact_suffix_trie(/* args */);
};

compact_suffix_trie::compact_suffix_trie(/* args */) {}