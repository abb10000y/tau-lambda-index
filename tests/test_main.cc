#include <gtest/gtest.h>
#include "symbol_table_test.h"
#include "k_factor_tree_test.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
