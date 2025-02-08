#define TESTING

#include <gtest/gtest.h>
#include "../include/Builder.h"
#include "boost/dynamic_bitset.hpp"
#include <random>
#include <algorithm>

class BuilderTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        sorted_keys = {1, 2, 3, 4, 5, 6, 7, 8};
        sparse_keys = {3, 7, 15, 31, 63, 127, 255, 511};
        sequential_keys = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
        normal_keys = {6, 8, 10, 11, 14, 16, 19, 29, 33, 34, 52, 57, 60, 62, 69,
                       73, 76, 81, 90, 92, 98, 110, 114, 117, 123, 126, 129, 150, 151, 156,
                       164, 166, 167, 170, 171, 172, 173, 174, 184, 188, 197, 199, 201, 212,
                       217, 219, 226, 228, 235, 242};
    }

    std::vector<uint32_t> sorted_keys;
    std::vector<uint32_t> sparse_keys;
    std::vector<uint32_t> sequential_keys;
    std::vector<uint64_t> normal_keys;
};

// Constructor Tests
TEST_F(BuilderTest, ConstructorValidation)
{
    EXPECT_NO_THROW(Builder<uint32_t>(sorted_keys, 4, 2));
    EXPECT_NO_THROW(Builder<uint32_t>(sparse_keys, 8, 4));
}

TEST_F(BuilderTest, ConstructorInitialization)
{
    Builder<uint32_t> builder(sorted_keys, 4, 2);

    EXPECT_EQ(builder.num_bins_, 4);
    EXPECT_EQ(builder.log_num_bins_, 2);
    EXPECT_EQ(builder.max_error_, 2);
    EXPECT_EQ(builder.num_keys_, 8);
    EXPECT_EQ(builder.min_key_, 1);
    EXPECT_EQ(builder.max_key_, 8);
    // 8 - 1 = 7 -> ceil log2(7) = 3 (log_range) -> 3 - 2 (log_num_bins) = 1 
    EXPECT_EQ(builder.shift_, 1);
    // 8 - 1 = 7 -> ceil log2(7) = 3 -> 2^3 = 8
    EXPECT_EQ(builder.range_, 8);
}

TEST_F(BuilderTest, ConstructorEdgeCases)
{
    // Invalid num_bins (not power of 2)
    EXPECT_DEATH(Builder<uint32_t>(sorted_keys, 3, 2), ".*");
    EXPECT_DEATH(Builder<uint32_t>(sorted_keys, 7, 2), ".*");

    // Invalid max_error (le 1)
    EXPECT_DEATH(Builder<uint32_t>(sorted_keys, 4, 1), ".*");
    EXPECT_DEATH(Builder<uint32_t>(sorted_keys, 4, 0), ".*");

    // Empty keys
    EXPECT_DEATH(Builder<uint32_t>({}, 4, 2), ".*");

    // num_bins > num_keys
    EXPECT_DEATH(Builder<uint32_t>({1, 2}, 4, 2), ".*");

    // max_error >= num_keys
    EXPECT_DEATH(Builder<uint32_t>({1, 2}, 4, 4), ".*");
}

// Bit Vector Creation Tests
TEST_F(BuilderTest, CreateBitVectorBasic)
{
    Builder<uint32_t> builder(sorted_keys, 4, 2);
    auto bit_vector = builder.createBitVector(sorted_keys);

    EXPECT_EQ(bit_vector.size(), 8);
    for (uint32_t key : sorted_keys)
    {
        EXPECT_TRUE(bit_vector[key - builder.min_key_]);
    }
}

TEST_F(BuilderTest, CreateBitVectorSparse)
{
    Builder<uint32_t> builder(sparse_keys, 4, 2);
    auto bit_vector = builder.createBitVector(sparse_keys);

    EXPECT_EQ(bit_vector.size(), 512);
    for (uint32_t key : sparse_keys)
    {
        EXPECT_TRUE(bit_vector[key - builder.min_key_]);
    }

    // Check random positions that should be empty
    std::vector<uint32_t> empty_positions = {10, 20, 100, 200, 300, 400};
    for (auto pos : empty_positions)
    {
        if (std::find(sparse_keys.begin(), sparse_keys.end(), pos) == sparse_keys.end())
        {
            EXPECT_FALSE(bit_vector[pos - builder.min_key_]);
        }
    }
}

TEST_F(BuilderTest, CreateBitVectorSIMDComparison)
{
    Builder<uint32_t> builder(sorted_keys, 4, 2);

    auto regular_vector = builder.createBitVector(sorted_keys);
    auto simd_vector = builder.createBitVectorSIMD(sorted_keys);

    EXPECT_EQ(regular_vector, simd_vector);
}

// Partition Tests
TEST_F(BuilderTest, PartitionVectorUniform)
{
    Builder<uint32_t> builder(sequential_keys, 4, 2);
    boost::dynamic_bitset<> bit_vector(16);
    for (size_t i = 0; i < 16; ++i)
    {
        bit_vector[i] = 1;
    }

    auto bins = builder.partitionVectorSIMD(bit_vector);

    EXPECT_EQ(bins.size(), 4);
    for (const auto &bin : bins)
    {
        EXPECT_EQ(bin.size(), 4);
        EXPECT_EQ(bin.count(), 4);
    }
}

TEST_F(BuilderTest, PartitionVectorUneven)
{
    Builder<uint32_t> builder(sorted_keys, 4, 2);
    auto bit_vector = builder.createBitVector(sorted_keys);
    auto bins = builder.partitionVectorSIMD(bit_vector);

    EXPECT_EQ(bins.size(), 4);
    size_t total_count = 0;
    for (const auto &bin : bins)
    {
        total_count += bin.count();
    }
    EXPECT_EQ(total_count, sorted_keys.size());
}

TEST_F(BuilderTest, PartitionVectorPrecise)
{
    // dynamic_bitset is little endian
    boost::dynamic_bitset<> bit_vector(std::string("01110101")); // {true, false, true, false, true, true, true, false}
    Builder<uint32_t> builder({1, 2, 3, 4}, 2, 3);

    std::vector<boost::dynamic_bitset<>> bins = builder.partitionVectorSIMD(bit_vector);

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

// Bin Counting Tests
TEST_F(BuilderTest, CountBinElementsUniform)
{
    Builder<uint32_t> builder(sequential_keys, 4, 2);
    std::vector<boost::dynamic_bitset<>> bins;
    for (int i = 0; i < 4; ++i)
    {
        boost::dynamic_bitset<> bin(8);
        for (int j = 0; j < 4; ++j)
        {
            bin[j] = 1;
        }
        bins.push_back(bin);
    }

    auto counts = builder.countBinElements(bins);

    EXPECT_EQ(counts.size(), 4);
    for (auto count : counts)
    {
        EXPECT_EQ(count, 4);
    }
}

TEST_F(BuilderTest, CountBinElementsEmpty)
{
    Builder<uint32_t> builder(sequential_keys, 4, 2);
    std::vector<boost::dynamic_bitset<>> bins(4, boost::dynamic_bitset<>(8));

    auto counts = builder.countBinElements(bins);

    EXPECT_EQ(counts.size(), 4);
    for (auto count : counts)
    {
        EXPECT_EQ(count, 0);
    }
}

// Build Tests
TEST_F(BuilderTest, BuildLeafOnlyTree)
{
    // Small dataset that should result in a leaf-only tree
    std::vector<uint32_t> small_keys = {1, 2, 3, 4};
    Builder<uint32_t> builder(small_keys, 4, 4);
    auto tree = builder.build();

    EXPECT_TRUE(tree.inner_nodes_.empty());
    EXPECT_FALSE(tree.leaf_nodes_.empty());
    EXPECT_EQ(tree.leaf_nodes_.size(), 4);
    for (size_t i = 0; i < 4; ++i)
    {
        EXPECT_EQ(tree.leaf_nodes_[i], 1);
    }
}

// Handchecked tree
TEST_F(BuilderTest, BuildNormalTree)
{
    Builder<uint64_t> builder(normal_keys, 4, 4);
    auto tree = builder.build();

    // Verify inner nodes
    std::vector<uint32_t> expected_inner_nodes = {
        15, 12, 14, 9, 8, 2147483648, 16, 24, 7, 3, 1, 4, 2147483652,
        4294967295, 4294967295, 2147483656, 0, 4, 7, 3, 4294967295,
        2147483660, 32, 4294967295, 3, 4, 2, 0, 4294967295, 2147483664,
        4294967295, 4294967295, 2, 4, 1, 0, 4294967295, 2147483668,
        4294967295, 4294967295};

    ASSERT_EQ(tree.inner_nodes_.size(), expected_inner_nodes.size());
    for (size_t i = 0; i < expected_inner_nodes.size(); ++i)
    {
        EXPECT_EQ(tree.inner_nodes_[i], expected_inner_nodes[i])
            << "Mismatch at inner node index " << i;
    }

    // Verify leaf nodes
    std::vector<uint32_t> expected_leaf_nodes = {
        3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 0, 1, 1, 1, 0, 2,
        1, 1, 1, 1};

    ASSERT_EQ(tree.leaf_nodes_.size(), expected_leaf_nodes.size());
    for (size_t i = 0; i < expected_leaf_nodes.size(); ++i)
    {
        EXPECT_EQ(tree.leaf_nodes_[i], expected_leaf_nodes[i])
            << "Mismatch at leaf node index " << i;
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}