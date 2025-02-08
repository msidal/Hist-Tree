#define TESTING

#include <gtest/gtest.h>
#include "../include/HistTree.h"
#include "../include/Builder.h"
#include <random>
#include <algorithm>

class HistTreeTest : public ::testing::Test
{
protected:
    HistTreeTest() : base_tree(createBaseTree()) {}

    void SetUp() override
    {
        normal_keys = {6, 8, 10, 11, 14, 16, 19, 29, 33, 34, 52, 57, 60, 62, 69,
                       73, 76, 81, 90, 92, 98, 110, 114, 117, 123, 126, 129, 150, 151, 156,
                       164, 166, 167, 170, 171, 172, 173, 174, 184, 188, 197, 199, 201, 212,
                       217, 219, 226, 228, 235, 242};
    }

    static HistTree<uint64_t> createBaseTree()
    {
        std::vector<uint64_t> init_keys = {6, 8, 10, 11, 14, 16, 19, 29, 33, 34, 52, 57, 60, 62, 69,
                                           73, 76, 81, 90, 92, 98, 110, 114, 117, 123, 126, 129, 150, 151, 156,
                                           164, 166, 167, 170, 171, 172, 173, 174, 184, 188, 197, 199, 201, 212,
                                           217, 219, 226, 228, 235, 242};
        Builder<uint64_t> builder(init_keys, 4, 4);
        return builder.build();
    }

    // Helper function to check if a key exists in the tree using bit_vector
    bool keyExistsInTree(const HistTree<uint64_t> &tree, uint64_t key)
    {
        if (key < tree.min_key_ || key > tree.max_key_)
            return false;
        return tree.bit_vector_[key - tree.min_key_];
    }

    std::vector<uint64_t> normal_keys;
    HistTree<uint64_t> base_tree;
};

// Insert Tests
TEST_F(HistTreeTest, InsertExistingKey)
{
    uint64_t existing_key = normal_keys[10];
    EXPECT_TRUE(base_tree.insert(existing_key));
    EXPECT_TRUE(keyExistsInTree(base_tree, existing_key));
}

TEST_F(HistTreeTest, InsertNewKeyInRange)
{
    uint64_t new_key = 15; // Between 14 and 16
    EXPECT_TRUE(base_tree.insert(new_key));
    EXPECT_TRUE(keyExistsInTree(base_tree, new_key));

    EXPECT_EQ(base_tree.num_keys_, normal_keys.size() + 1);
}

TEST_F(HistTreeTest, InsertNewKeyOutOfRange)
{
    uint64_t new_key = 300;
    EXPECT_TRUE(base_tree.insert(new_key));
    EXPECT_TRUE(keyExistsInTree(base_tree, new_key));

    EXPECT_EQ(base_tree.range_, 512);
}

// Remove Tests
TEST_F(HistTreeTest, RemoveExistingKey)
{
    uint64_t key_to_remove = normal_keys[20];
    EXPECT_TRUE(base_tree.remove(key_to_remove));
    EXPECT_FALSE(keyExistsInTree(base_tree, key_to_remove));

    EXPECT_EQ(base_tree.num_keys_, normal_keys.size() - 1);
}

TEST_F(HistTreeTest, RemoveNonExistentKey)
{
    uint64_t non_existent_key = 1000;
    EXPECT_FALSE(base_tree.remove(non_existent_key));
}

TEST_F(HistTreeTest, RemoveBoundaryKeys)
{
    uint64_t min_key = base_tree.min_key_;
    EXPECT_TRUE(base_tree.remove(min_key));
    EXPECT_FALSE(keyExistsInTree(base_tree, min_key));
    EXPECT_EQ(base_tree.min_key_, 8);

    uint64_t max_key = base_tree.max_key_;
    EXPECT_TRUE(base_tree.remove(max_key));
    EXPECT_FALSE(keyExistsInTree(base_tree, max_key));
    EXPECT_EQ(base_tree.max_key_, 235);
}

// GetSearchBound Tests
TEST_F(HistTreeTest, GetSearchBoundExistingKey)
{
    for (uint64_t key : normal_keys)
    {
        auto bound = base_tree.getSearchBound(key);
        EXPECT_TRUE(bound.end - bound.start <= base_tree.max_error_ + 1);

        // The key should exist in the tree and be within the bound
        bool found = false;
        for (size_t i = bound.start; i < bound.end; ++i)
        {
            if(normal_keys[i] == key)
            {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }
}

TEST_F(HistTreeTest, GetSearchBoundNonExistentKey)
{
    uint64_t non_existent_key = 15; // Between 14 and 16
    auto bound = base_tree.getSearchBound(non_existent_key);

    EXPECT_TRUE(bound.end - bound.start <= base_tree.max_error_ + 1);
    EXPECT_FALSE(keyExistsInTree(base_tree, non_existent_key));

    // Check if next greater key is within the bound
    bool found = false;
    for (size_t i = bound.start; i < bound.end; ++i)
    {
        if (normal_keys[i] == 16)
        {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
}

TEST_F(HistTreeTest, GetSearchBoundOutOfRange)
{
    // For key < min_key, lookup returns 0 and we get a normal bound
    uint64_t too_small = base_tree.min_key_ - 1;
    auto bound_small = base_tree.getSearchBound(too_small);
    EXPECT_EQ(bound_small.start, 0);
    EXPECT_GT(bound_small.end, 0);
    EXPECT_TRUE(bound_small.end <= base_tree.max_error_ + 1);

    // For key > max_key, lookup throws which is caught and returns {0,0}
    uint64_t too_large = base_tree.max_key_ + 10;
    auto bound_large = base_tree.getSearchBound(too_large);
    EXPECT_EQ(bound_large.start, 0);
    EXPECT_EQ(bound_large.end, 0);
}

// Stress Test
TEST_F(HistTreeTest, StressTest)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(1, 1000);
    std::uniform_int_distribution<int> op_dis(0, 1); // 0 for insert, 1 for remove

    Builder<uint64_t> stress_builder(normal_keys, 4, 4);
    HistTree<uint64_t> tree = stress_builder.build();

    // Perform 1000 random operations
    for (int i = 0; i < 1000; ++i)
    {
        uint64_t key = dis(gen);
        bool exists_before = keyExistsInTree(tree, key);

        if (op_dis(gen) == 0)
        { // Insert
            tree.insert(key);
            EXPECT_TRUE(keyExistsInTree(tree, key));
        }
        else
        { // Remove
            bool remove_success = tree.remove(key);
            EXPECT_EQ(remove_success, exists_before);
            EXPECT_FALSE(keyExistsInTree(tree, key));
        }
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}