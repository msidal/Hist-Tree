#define TESTING

#include <gtest/gtest.h>
#include "../include/HistTree.h"

TEST(HistTreeTest, dummy) {
    ASSERT_EQ(1, 1);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}