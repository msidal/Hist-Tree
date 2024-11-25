#include "Builder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

#include <HistTree.h>

template <typename KeyType>
Builder<KeyType>::Builder(KeyType min_key, KeyType max_key, size_t num_bins, size_t max_error) 
    :   min_key_(min_key),
        max_key_(max_key), 
        num_bins_(num_bins), 
        log_num_bins_(computeLog(static_cast<uint64_t>(num_bins))),
        max_error_(max_error),
        curr_num_keys_(0),
        prev_key_(min_key) {
    assert((num_bins & (num_bins - 1)) == 0); // num_bins must be a power of 2
    assert(min_key < max_key);
    assert(num_bins > 1);
    assert(max_error > 0);
    assert(max_error < num_bins); 

    // logarithm of the range of keys 
    auto lg = computeLog(max_key - min_key, true); 

    assert(lg >= log_num_bins_);
    shift_ = lg - log_num_bins_;
}

template <typename KeyType>
void Builder<KeyType>::insert(KeyType key) {
    assert(key >= min_key_);
    assert(key <= max_key_);
    assert(key >= prev_key_);

    curr_num_keys_++;
    


    prev_key_ = key;
}
