#pragma once

#include <vector>
#include <immintrin.h>
#include <cassert>
#include <cmath>
#include <limits>

#include <common.h>
#include <HistTree.h>


template <typename KeyType>
class Builder {
    public:
        Builder(KeyType min_key, KeyType max_key, size_t num_bins, size_t max_error, vector<KeyType> keys) 
            :   min_key_(min_key),
                max_key_(max_key), 
                num_bins_(num_bins), 
                log_num_bins_(computeLog(static_cast<uint64_t>(num_bins))),
                max_error_(max_error),
                curr_num_keys_(0),
                prev_key_(min_key),
                keys_(std::move(keys)) {
            assert((num_bins & (num_bins - 1)) == 0); // num_bins must be a power of 2
            assert(min_key < max_key);
            assert(num_bins > 1);
            assert(max_error > 0);
            assert(max_error < num_bins); 

            max_height_ = computeMaxHeight();

            // for efficient computation
            shift_ = computeLog(num_bins);

            // range_ rounds up to the next power of 2 
            range_ = 1 << static_cast<unsigned>(std::ceil(std::log2(max_key - min_key)));
        } 

        HistTree<KeyType> build() {
            auto bit_vector = createBitVector(keys_);
            auto bins = partitionBitVector(bit_vector, num_bins_);

            // count the number of set bits in each partition
            std::vector<size_t> counts(num_bins_, 0);
            for (size_t i = 0; i < num_bins_; ++i) {
                for (uint32_t value : partitions[i]) {
                    partition_counts[i] += countSetBits(value);
                }   
            }

            buildNodes(partition_counts, partitions);

            return HistTree<KeyType>(min_key_, max_key_, keys_.size(), num_bins_, log_num_bins_, max_error_, shift_, inner_nodes_, leaf_nodes_);
        }

    private:
        void buildNodes(const std::vector<size_t>& partition_counts,
                const std::vector<std::vector<uint32_t>>& partitions) {
            for (size_t i = 0; i < partition_counts.size(); ++i) {
                int count = partition_counts[i];
                if (count >= max_error_) {
                    // Wenn die Zählung den Threshold überschreitet, erstelle einen neuen inneren Knoten
                    inner_nodes_.push_back(count);

                    // Partitioniere die aktuelle Partition erneut, da sie den Threshold überschreitet
                    std::vector<std::vector<uint32_t>> new_partitions = partitionBins(partitions[i], num_bins_);
                    std::vector<size_t> new_partition_counts(num_bins_, 0);

                    // Zähle die Bits in den neuen Partitionen
                    for (size_t j = 0; j < num_bins_; ++j) {
                        for (uint32_t value : new_partitions[j]) {
                            new_partition_counts[j] += countSetBits(value);
                        }
                    }

                    // Rekursiv weiter partitionieren, wenn nötig
                    buildNodes(new_partition_counts, new_partitions);
                } else {
                    // Wenn die Zählung den Threshold nicht überschreitet, erstelle einen Blattknoten
                    leaf_nodes_.push_back(count);
                }
            }
        }

        // create a bit vector from the keys
        std::vector<uint32_t> createBitVector(const std::vector<KeyType>& keys) {
            size_t bit_length = (range_ + 31) / 32;
            std::vector<uint32_t> bit_vector(bit_length, 0);

            for (auto key : keys) {
                int offset = key - min_key_;
                bit_vector[offset / 32] |= (1 << (offset % 32));
            }

            return bit_vector;
        }

        // Only works for 32-bit integers 
        // g++ -O2 -mavx2 -mfma -std=c++17 -o simd_example main.cpp
        std::vector<uint32_t> createBitVectorSIMD(const std::vector<uint32_t>& keys) {
            size_t bit_length = (range_ + 31) / 32; // amount of 32-bit blocks
            std::vector<uint32_t> bit_vector(bit_length, 0);

            // SIMD-specific parameters
            const size_t SIMD_WIDTH = 8; // amount of parallel processing (256-bit / 32-bit = 8)
            __m256i vec_min_value = _mm256_set1_epi32(min_value); // SIMD-vector for min_value

            size_t i = 0;
            for (; i + SIMD_WIDTH <= keys.size(); i += SIMD_WIDTH) {
                // load 8 keys into a SIMD register
                __m256i vec_keys = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&keys[i]));

                // subtract min_value from each key
                __m256i vec_offsets = _mm256_sub_epi32(vec_keys, vec_min_value);

                // block_index = offsets / 32
                __m256i vec_block_index = _mm256_srli_epi32(vec_offsets, 5); // Division durch 32 durch Rechts-Shift

                // bit_offset = offsets % 32
                __m256i vec_bit_offset = _mm256_and_si256(vec_offsets, _mm256_set1_epi32(31)); // Modulo 32 mit Maske

                // save SIMD registers in arrays
                alignas(32) int block_indices[SIMD_WIDTH];
                alignas(32) int bit_offsets[SIMD_WIDTH];

                _mm256_storeu_si256(reinterpret_cast<__m256i*>(block_indices), vec_block_index);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(bit_offsets), vec_bit_offset);

                // Setze die Bits im Bitvektor
                for (int j = 0; j < SIMD_WIDTH; ++j) {
                    if (block_indices[j] < static_cast<int>(bit_vector.size())) {
                        bit_vector[block_indices[j]] |= (1 << bit_offsets[j]);
                    }
                }
            }

            // remaining keys (if keys.size() is not divisible by SIMD_WIDTH)
            for (; i < keys.size(); ++i) {
                int index = keys[i] - min_key_;
                if (index >= 0 && static_cast<size_t>(index) < range_) {
                    bit_vector[index / 32] |= (1 << (index % 32));
                }
            }

            return bit_vector;
        }

        // partition the bit vector into bins
        std::vector<std::vector<uint32_t>> partitionBitVector(const std::vector<uint32_t>& bit_vector, size_t num_bins) {
            // calculate the total number of bits
            size_t total_bits = bit_vector.size() * 32;

            // calculate how many bits per bin
            size_t bits_per_bin = (total_bits + num_bins - 1) / num_bins;

            std::vector<std::vector<uint32_t>> bins(num_bins);

            size_t current_bit = 0;
            size_t bin_index = 0;

            // go through each block in the bit vector
            for (size_t block_index = 0; block_index < bit_vector.size(); ++block_index) {
                uint32_t block = bit_vector[block_index];
                
                // determine the bin for this block
                while (current_bit >= (bin_index + 1) * bits_per_bin) {
                    ++bin_index;
                }
                
                // add the block to the bin
                bins[bin_index].push_back(block);
                
                current_bit += 32; // move to the next block
            }

            return bins;
        }

        static int countSetBits(uint32_t) {
            return __builtin_popcount(n);
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
        constexpr int32_t setHighOrderBit(int32_t value) {
            return value | (1 << 31);
        }

        constexpr bool isHighOrderBitSet(int32_t value) {
            return value & (1 << 31);
        }

        constexpr int32_t clearHighOrderBit(int32_t value) {
            return value & ~(1 << 31);
        }

        static constexpr unsigned Terminal = 0xFFFFFFFF; // marks that the bin is terminal

        const KeyType min_key_;
        const KeyType max_key_;
        const size_t num_bins_;
        const size_t log_num_bins_;
        const size_t max_error_;
        //const bool single_pass_;
        //const bool use_cache_;
        const size_t range_;
        const size_t max_height_;

        size_t curr_num_keys_;
        KeyType prev_key_;
        size_t shift_;

        std::vector<KeyType> keys_;
        std::vector<uint32_t> inner_nodes_;
        std::vector<uint32_t> leaf_nodes_;

};