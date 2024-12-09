#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <memory>
#include <iostream>

#include <common.h>

template <typename KeyType>
class HistTree {
public:
    HistTree() = default; 

    HistTree(KeyType min_key, KeyType max_key, size_t num_keys, size_t num_bins, 
            size_t log_num_bins, size_t max_error, size_t shift, size_t range,
            std::vector<uint32_t> inner_nodes, std::vector<uint32_t> leaf_nodes) :
            min_key_(min_key),
            max_key_(max_key),
            num_keys_(num_keys),
            num_bins_(num_bins),
            log_num_bins_(log_num_bins),
            max_error_(max_error),
            shift_(shift),
            range_(range),
            inner_nodes_(std::move(inner_nodes)),
            leaf_nodes_(std::move(leaf_nodes)) {};

    SearchBound getSearchBound( const KeyType key) const {
        size_t begin = 0;
        try {
            begin = lookup(key);
        } catch (const std::out_of_range& e) {
            std::cerr << e.what() << std::endl;
            return SearchBound{0, 0};
        }
        // `end` is exclusive.
        const size_t end = (begin + max_error_ + 1 > num_keys_)
                            ? num_keys_
                            : (begin + max_error_ + 1);
        return SearchBound{begin, end};
    }

    // lookup of first key equal or greater than key
    size_t lookup(KeyType key) const {
        if (key < min_key_) return 0;
        if (key > max_key_) throw std::out_of_range("key out of range");
        
        // if root is a leaf node
        if (inner_nodes_.empty()) {
            
        }

        uint32_t next, *node = const_cast<uint32_t*>(inner_nodes_.data());
        bool done = false;
        size_t width = shift_;
        size_t i, bin, pos = 0;
        key -= min_key_;

        do {
            bin = key >> width;
            for(i = 0; i < bin; i++) {
                pos += node[i];
            }

            if (done || node[num_bins_ + bin] == Terminal) return pos;

            done = isHighOrderBitSet(node[num_bins_ + bin]);
            next = clearHighOrderBit(node[num_bins_ + bin]);
            node = done ? const_cast<uint32_t*>(leaf_nodes_.data()) + next 
                        : const_cast<uint32_t*>(inner_nodes_.data()) + next;

            key -= bin << width;
            width -= log_num_bins_;

        } while (1);
    }

    void insert(KeyType key) {

    }

    void remove(KeyType key) {

    }

    bool isEmpty() const {
        return num_keys_ == 0;
    }

    void print() const {
        
    };

    size_t getSize() const {
        return sizeof(*this) + inner_nodes_.size() * sizeof(uint32_t) + leaf_nodes_.size() * sizeof(uint32_t);
    }

#ifdef TESTING
    public:  
#else
    private: 
#endif
    static constexpr unsigned Terminal = 0xFFFFFFFF; // marks that the bin is terminal

    // helper functions for manipulating the high order bit
    constexpr uint32_t setHighOrderBit(uint32_t value) {
        return value | (1 << 31);
    }

    constexpr bool isHighOrderBitSet(uint32_t value) const {
        return value & (1 << 31);
    }

    constexpr uint32_t clearHighOrderBit(uint32_t value) const{
        return value & ~(1 << 31);
    }

     // helper functions for computing the log base 2 of a number
    static unsigned computeLog(uint32_t n, bool round = false) {
        assert(n);
        return 31 - __builtin_clz(n) + (round ? ((n & (n - 1)) != 0) : 0);
    }

    static unsigned computeLog(uint64_t n, bool round = false) {
        assert(n);
        return 63 - __builtin_clzl(n) + (round ? ((n & (n - 1)) != 0) : 0);
    }



    KeyType min_key_;
    KeyType max_key_;
    size_t num_keys_;
    const size_t range_;
    const size_t num_bins_;
    const size_t log_num_bins_;
    const size_t max_error_;
    const size_t shift_;

    // physical layout of the tree
    std::vector<uint32_t> inner_nodes_;
    std::vector<uint32_t> leaf_nodes_;
};
