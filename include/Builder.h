#pragma once

#include <cassert>
#include <cmath>
#include <limits>

#include <common.h>
#include <HistTree.h>


template <typename KeyType>
class Builder {
    public:
        Builder(KeyType min_key, KeyType max_key, size_t num_bins, size_t max_error); // TODO: add single_pass and use_cache

        void insert(KeyType key);

        HistTree<KeyType> build(); 
    private:
        static constexpr unsigned Infinity = std::numeric_limits<unsigned>::max();
        static constexpr unsigned Leaf = 1 << 31;

        // Range covered by a node, i.e. [l, r[
        using Range = std::pair<unsigned, unsigned>;

        // (Node level, smallest key in node)
        using Info = std::pair<unsigned, KeyType>;

        // A queue element
        using Elem = std::pair<unsigned, Range>;

         // helper functions for computing the log 
        static unsigned computeLog(uint32_t n, bool round = false) {
            assert(n);
            return 31 - __builtin_clz(n) + (round ? ((n & (n - 1)) != 0) : 0);
        }

        static unsigned computeLog(uint64_t n, bool round = false) {
            assert(n);
            return 63 - __builtin_clzl(n) + (round ? ((n & (n - 1)) != 0) : 0);
        }

        // helper functions for manipulating the high order bit
        constexpr int32_t setHighOrderBit(int32_t value) {
            return value | (1 << 31);
        }

        constexpr bool isHighOrderBitSet(int32_t value) {
            return value & (1 << 31);
        }

        constexpr int32_t clearHighOrderBit(int32_t value) {
            return value & ~(1 << 31);
        }

        void prepareTree();
        void prune();

        const KeyType min_key_;
        const KeyType max_key_;
        const size_t num_bins_;
        const size_t log_num_bins_;
        const size_t max_error_;
        const bool single_pass_;
        const bool use_cache_;

        size_t curr_num_keys_;
        KeyType prev_key_;
        size_t shift_;

        std::vector<KeyType> keys_;
        std::vector<uint32_t> inner_nodes_;
        std::vector<uint32_t> leaf_nodes_;
        std::vector<std::pair<Info, std::vector<Range>>> tree_;
};