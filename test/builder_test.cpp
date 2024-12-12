#define TESTING

#include <gtest/gtest.h>
#include "../include/Builder.h"

// subject to extension


TEST(BuilderTest, TestConstructor) {
    std::vector<uint32_t> keys = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    Builder<uint32_t>* builder = new Builder<uint32_t>(keys, 4, 2);

    ASSERT_EQ(builder->num_bins_, 4);
    ASSERT_EQ(builder->log_num_bins_, 2);
    ASSERT_EQ(builder->max_error_, 2);
    ASSERT_EQ(builder->num_keys_, 10);
    ASSERT_EQ(builder->min_key_, 1);
    ASSERT_EQ(builder->max_key_, 10);
    ASSERT_EQ(builder->shift_, 2);
    ASSERT_EQ(builder->range_, 16);
}

TEST(BuilderTest, TestFalseConstructor) {
    ASSERT_DEATH(new Builder<uint32_t>({1, 2, 3, 4}, 3, 2), ".*");
}

TEST(BuilderTest, ComputeLog) {
    EXPECT_EQ(Builder<uint32_t>::computeLog(static_cast<uint32_t>(2)), 1);
    EXPECT_EQ(Builder<uint32_t>::computeLog(static_cast<uint32_t>(4)), 2);
    EXPECT_EQ(Builder<uint32_t>::computeLog(static_cast<uint32_t>(8)), 3);
}


// high order bit manipulation cannot be tested as the functions are constexpr

TEST(BuilderTest, CreateBitVector) {
    std::vector<uint32_t> keys = {3, 7, 9, 10, 15, 26, 28, 30, 31, 40};
    Builder<uint32_t> builder(keys, 2, 2);

    std::vector<bool> bit_vector = builder.createBitVector(keys);

    ASSERT_EQ(bit_vector.size(), 64);
    
    // check if all keys are set
    for (uint32_t key : keys) {
        ASSERT_EQ(bit_vector[key - builder.min_key_], true);
    }
    // check if rest is not set
    for (uint32_t i = builder.min_key_; i <= builder.max_key_; ++i) {
        if (std::find(keys.begin(), keys.end(), i) == keys.end()) {
            ASSERT_EQ(bit_vector[i - builder.min_key_], false);
        }
    }
}

TEST(BuilderTest, CountBinElements) {
    std::vector<std::vector<bool>> bins = {
        {true, false, true, false},
        {true, true, true, false},
        {false, false, false, false},
        {true, true, true, true}
    };

    Builder<uint32_t> builder({1, 2, 3, 4}, 4, 2);

    std::vector<uint32_t> counts = builder.countBinElements(bins);

    ASSERT_EQ(counts.size(), 4);
    ASSERT_EQ(counts[0], 2);
    ASSERT_EQ(counts[1], 3);
    ASSERT_EQ(counts[2], 0);
    ASSERT_EQ(counts[3], 4);
}

TEST(BuilderTest, PartitionVector) {
    std::vector<bool> bit_vector = {true, false, true, false, true, true, true, false};
    Builder<uint32_t> builder({1, 2, 3, 4}, 2, 2);

    std::vector<std::vector<bool>> bins = builder.partitionVector(bit_vector);

    ASSERT_EQ(bins.size(), 2);
    ASSERT_EQ(bins[0].size(), 4);
    ASSERT_EQ(bins[1].size(), 4);
    ASSERT_EQ(bins[0][0], true);
    ASSERT_EQ(bins[0][1], false);
    ASSERT_EQ(bins[0][2], true);
    ASSERT_EQ(bins[0][3], false);
    ASSERT_EQ(bins[1][0], true);
    ASSERT_EQ(bins[1][1], true);
    ASSERT_EQ(bins[1][2], true);
    ASSERT_EQ(bins[1][3], false);
}

// TODO test build & createBitVectorSIMD

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}