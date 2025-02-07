#define TESTING

#include <gtest/gtest.h>
#include "../include/Builder.h"
#include "../include/HistTree.h"
#include <vector>
#include <algorithm>

//reset
TEST(HistTreeTest, DUMMY) {
    EXPECT_EQ(1, 1);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}