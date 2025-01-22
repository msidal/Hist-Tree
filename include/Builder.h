#pragma once

#include <vector>
#include <immintrin.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <queue>
#include <boost/dynamic_bitset.hpp>

#include "HistTree.h"


template <typename KeyType>
class Builder {
    public:
        Builder(std::vector<KeyType> keys, size_t num_bins, size_t max_error)
            : num_bins_(num_bins),
            log_num_bins_(computeLog(static_cast<uint64_t>(num_bins))),
            max_error_(max_error),
            num_keys_(keys.size()),
            keys_(std::move(keys)) {

            assert(keys_.size() > 0 && "keys should not be empty");
            assert(num_bins_ <= num_keys_ && "num_bins should not be greater than num_keys");

            min_key_ = keys_.front();
            max_key_ = keys_.back();
            assert(min_key_ < max_key_);
            assert((num_bins_ & (num_bins_ - 1)) == 0 && "num_bins should be a power of 2"); 
            assert(num_bins_ > 1 && "num_bins should be greater than 1");
            assert(max_error_ > 1 && "max_error should be greater than 1");

            // calculate the range and shift
            auto log_range_ = computeLog(max_key_ - min_key_, true);
            assert(log_range_ >= log_num_bins_ && "range should be greater than log_num_bins");
            shift_ = log_range_ - log_num_bins_;
            range_ = 1 << log_range_;
            
        }

        HistTree<KeyType> build() {
            // estimated size of the tree
            size_t estimated_max_depth = static_cast<size_t>(std::ceil(std::log(static_cast<double>(num_keys_) / max_error_) / std::log(num_bins_)));
            size_t estimated_inner_size = ((num_keys_ / max_error_) * (estimated_max_depth + 1) - 1) << (log_num_bins_ + 1);
            size_t estimated_leaf_size = static_cast<size_t>(std::pow(num_bins_, estimated_max_depth)) << log_num_bins_;

            inner_nodes_.resize(estimated_inner_size);
            leaf_nodes_.resize(estimated_leaf_size);

            //auto bit_vector = createBitVector(keys_);
            auto bit_vector = createBitVectorSIMD(keys_);
            auto bins = partitionVector(bit_vector); 
            
            std::queue<std::tuple<std::vector<boost::dynamic_bitset<>>, size_t, size_t>> to_process;
            to_process.push({bins, 0, 0});

            size_t next_inner_index = 0;
            size_t next_leaf_index = 0;

            while(!to_process.empty()) {
                auto [bins, current_index, bin_index] = to_process.front();
                to_process.pop();

                // calculate the number of elements in each bin
                auto counts = countBinElements(bins);
                std::vector<uint32_t> child_count;

                // check if node is a leaf
                if (*std::max_element(counts.begin(), counts.end()) < max_error_) {
                    if (current_index > 0) {
                        inner_nodes_[current_index + num_bins_ + bin_index] = setHighOrderBit(next_leaf_index);

                    }
                    for (size_t i = 0; i < num_bins_; ++i) {
                        leaf_nodes_[next_leaf_index + i] = counts[i];
                    }
                    next_leaf_index += num_bins_;
                } else {
                    for (size_t i = 0; i < num_bins_; ++i) {
                        inner_nodes_[current_index + i] = counts[i];
                    }
                    // check if the current node is the root and set the next next free inner index accordingly
                    if (current_index == 0) next_inner_index = num_bins_ << 1;
                    // go through each bin and check if it is terminal or a new node by checking the count
                    for (size_t i = 0; i < num_bins_; ++i) {
                        // terminal bin
                        if (counts[i] < max_error_) {
                            inner_nodes_[current_index + num_bins_ + i] = Terminal;
                        } else {
                            auto child_bins = partitionVector(bins[i]);
                            child_count = countBinElements(child_bins);
                            //check if the child node is a leaf
                            if (*std::max_element(child_count.begin(), child_count.end()) < max_error_) {
                                inner_nodes_[current_index + num_bins_ + i] = setHighOrderBit(next_leaf_index);
                                for (size_t i = 0; i < counts.size(); ++i) {
                                    leaf_nodes_[next_leaf_index + i] = child_count[i];
                                }
                                next_leaf_index += num_bins_;
                            } else {
                                inner_nodes_[current_index + num_bins_ + i] = next_inner_index;
                                to_process.push({child_bins, next_inner_index, i});
                                next_inner_index += num_bins_ << 1;
                            }
                        }
                    }

                }
            }

            inner_nodes_.resize(next_inner_index); 
            leaf_nodes_.resize(next_leaf_index);   

            inner_nodes_.shrink_to_fit();      
            leaf_nodes_.shrink_to_fit();

            return HistTree<KeyType>(min_key_, max_key_, num_keys_, num_bins_, log_num_bins_, max_error_, shift_, range_, inner_nodes_, leaf_nodes_, bit_vector);
        }


    #ifdef TESTING
        public:  
    #else
        private: 
    #endif
        // count set bits in all bins
        std::vector<uint32_t> countBinElements(const std::vector<boost::dynamic_bitset<>>& bins) {
            std::vector<uint32_t> counts(num_bins_);

            for (size_t i = 0; i < num_bins_; ++i) {
                counts[i] = bins[i].count();
            }

            return counts;
        }

        // partition the bit vector into bins
        std::vector<boost::dynamic_bitset<>> partitionVector(const boost::dynamic_bitset<>& bitset) {
            size_t bin_size = bitset.size() / num_bins_;
            std::vector<boost::dynamic_bitset<>> result;
            size_t start_index = 0;

            for (size_t i = 0; i < num_bins_; ++i) {
                boost::dynamic_bitset<> bin(bin_size);

                for (size_t j = 0; j < bin_size; ++j) {
                    bin[j] = bitset[start_index + j];
                }

                result.push_back(bin);

                start_index += bin_size;
            }

            return result;
        }
        
        // create a bit vector
        boost::dynamic_bitset<> createBitVector(const std::vector<KeyType>& keys) {
            boost::dynamic_bitset<> bit_vector(range_);

            for (auto key : keys) {
                int offset = key - min_key_; // map key to an offset in the bit vector
                bit_vector.set(offset); 
            }

            return bit_vector;
        }

        // SIMD-optimized access for uint32_t keys
        boost::dynamic_bitset<> createBitVectorSIMD(const std::vector<KeyType>& keys) {
            boost::dynamic_bitset<> bit_vector(range_);
            __m256i min_key_vec = _mm256_set1_epi32(min_key_); // load min_key into a SIMD register
            
            for (size_t i = 0; i < keys.size(); i += 8) { // SIMD batch of 8 keys
                size_t remaining = keys.size() - i;
                int temp_keys[8] = {0};

                for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j) {
                    temp_keys[j] = keys[i + j];
                }

                __m256i key_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(temp_keys));
                __m256i offset_vec = _mm256_sub_epi32(key_vec, min_key_vec); // offset calculation

                alignas(32) int offsets[8];
                _mm256_store_si256(reinterpret_cast<__m256i*>(offsets), offset_vec); // save values

                for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j) {
                    if (offsets[j] >= 0 && offsets[j] < static_cast<int>(range_)) {
                        bit_vector.set(offsets[j]); // set bit
                    }
                }
            }
            return bit_vector;
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

        // helper functions for manipulating the high order bit
        constexpr uint32_t setHighOrderBit(uint32_t value) {
            return value | (1 << 31);
        }

        constexpr bool isHighOrderBitSet(uint32_t value) {
            return value & (1 << 31);
        }

        constexpr uint32_t clearHighOrderBit(uint32_t value) {
            return value & ~(1 << 31);
        }

        static constexpr unsigned Terminal = 0xFFFFFFFF; // marks that the bin is terminal

        KeyType min_key_;
        KeyType max_key_;
        const size_t num_bins_;
        const size_t log_num_bins_;
        const size_t max_error_;
        //const bool single_pass_;
        //const bool use_cache_;
        size_t range_;

        // for building the tree
        size_t next_inner_index = 0;
        size_t num_keys_;
        KeyType prev_key_;
        size_t shift_;

        std::vector<KeyType> keys_;
        std::vector<uint32_t> inner_nodes_;
        std::vector<uint32_t> leaf_nodes_;

};

