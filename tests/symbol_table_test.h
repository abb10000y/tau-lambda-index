#pragma once
#include <gtest/gtest.h>
#include "../include/symbol_table/symbol_table.h"

// 測試的函數
int Add(int a, int b) {
    return a + b;
}

// 測試案例
TEST(AdditionTest, PositiveNumbers) {
    EXPECT_EQ(Add(2, 3), 5);
}

TEST(AdditionTest, NegativeNumbers) {
    EXPECT_EQ(Add(-2, -3), -5);
}