#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <queue>
#include <memory>

#include <common.h>

template <typename KeyType>
class HistTree {
public:
    HistTree() = default; 

    HistTree(KeyType min_key, KeyType max_key, size_t num_keys, size_t num_bins, 
            size_t log_num_bins, size_t max_error, size_t shift, 
            std::vector<unsigned> inner_nodes, std::vector<unsigned> leaf_nodes);

    SearchBound getSearchBound( const KeyType key) const;

    size_t lookup(KeyType key) const;

    void insert(KeyType key, size_t value);

    void remove(KeyType key, size_t value);

    void update(KeyType key, size_t old_value, size_t new_value);

    bool isEmpty() const {
        return num_keys_ == 0;
    }

    void print() const;

    size_t getSize() const {
        return sizeof(*this) + inner_nodes.size() * sizeof(unsigned) + leaf_nodes.size() * sizeof(unsigned);
    }

private:
    KeyType min_key_;
    KeyType max_key_;
    size_t num_keys_;
    size_t num_bins_;
    size_t log_num_bins_;
    size_t max_error_;
    size_t shift_;

    //physical layout of the tree
    std::vector<uint32_t> inner_nodes_;
    std::vector<uint32_t> leaf_nodes_;
};
