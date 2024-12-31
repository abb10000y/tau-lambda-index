#include <gtest/gtest.h>

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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
