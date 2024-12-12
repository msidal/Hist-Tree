#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <queue>
#include <memory>
#include <iostream>

#include <common.h>

template <typename KeyType>
class HistTree {
public:
    HistTree() = default; 

    HistTree(KeyType min_key, KeyType max_key, size_t num_keys, size_t num_bins, 
            size_t log_num_bins, size_t max_error, size_t shift, 
            std::vector<uint32_t> inner_nodes, std::vector<uint32_t> leaf_nodes) :
            min_key_(min_key),
            max_key_(max_key),
            num_keys_(num_keys),
            num_bins_(num_bins),
            log_num_bins_(log_num_bins),
            max_error_(max_error),
            shift_(shift),
            inner_nodes_(std::move(inner_nodes)),
            leaf_nodes_(std::move(leaf_nodes)) {};

    SearchBound getSearchBound( const KeyType key) const {
        return {0, 0};
    }

    size_t lookup(KeyType key) const {
        return 0;
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

private: 
    static constexpr unsigned Leaf = (1u << 31);
    static constexpr unsigned Mask = Leaf - 1;

    KeyType min_key_;
    KeyType max_key_;
    size_t num_keys_;
    const size_t num_bins_;
    const size_t log_num_bins_;
    const size_t max_error_;
    const size_t shift_;

    //physical layout of the tree
    std::vector<uint32_t> inner_nodes_;
    std::vector<uint32_t> leaf_nodes_;
};
