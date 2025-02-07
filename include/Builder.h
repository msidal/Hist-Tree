#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <queue>
#include <boost/dynamic_bitset.hpp>
#include "omp.h"
#include "HistTree.h"

#ifdef __AVX2__
    #include <immintrin.h>
#elif defined(__ARM_NEON)
    #include <arm_neon.h>
#endif

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
            assert(max_error_ <= num_keys_ && "max_error should not be greater than num_keys");

            min_key_ = keys_.front();
            max_key_ = keys_.back();
            assert(min_key_ < max_key_);
            assert((num_bins_ & (num_bins_ - 1)) == 0 && "num_bins should be a power of 2"); 
            assert(num_bins_ > 1 && "num_bins should be greater than 1");
            assert(max_error_ > 1 && "max_error should be greater than 1");

            auto log_range_ = computeLog(max_key_ - min_key_, true);
            assert(log_range_ >= log_num_bins_ && "range should be greater than log_num_bins");
            shift_ = log_range_ - log_num_bins_;
            range_ = 1 << log_range_;
        }

        HistTree<KeyType> build() {
            size_t estimated_max_depth = static_cast<size_t>(std::ceil(std::log(static_cast<double>(num_keys_) / max_error_) / std::log(num_bins_)));
            size_t estimated_inner_size = ((num_keys_ / max_error_) * (estimated_max_depth + 1) - 1) << (log_num_bins_ + 1);
            size_t estimated_leaf_size = static_cast<size_t>(std::pow(num_bins_, estimated_max_depth)) << log_num_bins_;

            inner_nodes_.resize(estimated_inner_size);
            leaf_nodes_.resize(estimated_leaf_size);

            auto bit_vector = createBitVectorSIMD(keys_);
            auto bins = partitionVectorSIMD(bit_vector); 
            
            std::queue<std::tuple<std::vector<boost::dynamic_bitset<>>, size_t, size_t>> to_process;
            to_process.push({bins, 0, 0});

            size_t next_inner_index = 0;
            size_t next_leaf_index = 0;

            while(!to_process.empty()) {
                auto [bins, current_index, bin_index] = to_process.front();
                to_process.pop();

                auto counts = countBinElements(bins);
                std::vector<uint32_t> child_count;

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
                    if (current_index == 0) next_inner_index = num_bins_ << 1;
                    for (size_t i = 0; i < num_bins_; ++i) {
                        if (counts[i] < max_error_) {
                            inner_nodes_[current_index + num_bins_ + i] = Terminal;
                        } else {
                            auto child_bins = partitionVectorSIMD(bins[i]);
                            child_count = countBinElements(child_bins);
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

            return HistTree<KeyType>(min_key_, max_key_, num_keys_, num_bins_, log_num_bins_, max_error_, 
                                   shift_, range_, inner_nodes_, leaf_nodes_, bit_vector);
        }

    #ifdef TESTING
        public:  
    #else
        private: 
    #endif
        std::vector<uint32_t> countBinElements(const std::vector<boost::dynamic_bitset<>>& bins) {
            std::vector<uint32_t> counts(num_bins_);
            for (size_t i = 0; i < num_bins_; ++i) {
                counts[i] = bins[i].count();
            }
            return counts;
        }

        std::vector<boost::dynamic_bitset<>> partitionVector(const boost::dynamic_bitset<>& bitset) {
            size_t bin_size = bitset.size() / num_bins_;
            std::vector<boost::dynamic_bitset<>> result;
            result.reserve(num_bins_);
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

        std::vector<boost::dynamic_bitset<>> partitionVectorSIMD(const boost::dynamic_bitset<>& bitset) {
            const size_t bin_size = bitset.size() / num_bins_;
            std::vector<boost::dynamic_bitset<>> bins(num_bins_, boost::dynamic_bitset<>(bin_size));

            #if defined(__AVX2__) || defined(__ARM_NEON)
                #ifdef __AVX2__
                    constexpr size_t simd_width = 256;
                #else
                    constexpr size_t simd_width = 128;
                #endif
                constexpr size_t simd_bytes = simd_width / 8;
                const size_t L1_cache_size = 32768;
                const size_t chunk_size = (L1_cache_size / simd_bytes) * simd_width;

                #pragma omp parallel
                {
                    alignas(32) unsigned char simd_buffer[simd_bytes];
                    
                    #pragma omp for schedule(dynamic)
                    for (size_t bin = 0; bin < num_bins_; ++bin) {
                        const size_t offset = bin * bin_size;
                        
                        for (size_t chunk_start = 0; chunk_start < bin_size; chunk_start += chunk_size) {
                            const size_t current_chunk = std::min(chunk_size, bin_size - chunk_start);
                            
                            for (size_t i = 0; i < current_chunk; i += simd_width) {
                                const size_t bits_to_process = std::min(simd_width, current_chunk - i);
                                const size_t pos = offset + chunk_start + i;
                                
                                #ifdef __AVX2__
                                    __m256i bits = _mm256_setzero_si256();
                                    for (size_t j = 0; j < bits_to_process; j += 8) {
                                        unsigned char byte = 0;
                                        for (size_t k = 0; k < 8 && (j + k) < bits_to_process; ++k) {
                                            byte |= (bitset[pos + j + k] ? 1 : 0) << k;
                                        }
                                        simd_buffer[j/8] = byte;
                                    }
                                    bits = _mm256_load_si256(reinterpret_cast<const __m256i*>(simd_buffer));
                                    _mm256_store_si256(reinterpret_cast<__m256i*>(simd_buffer), bits);
                                #elif defined(__ARM_NEON)
                                    uint8x16_t bits = vdupq_n_u8(0);
                                    for (size_t j = 0; j < bits_to_process; j += 8) {
                                        uint8_t byte = 0;
                                        for (size_t k = 0; k < 8 && (j + k) < bits_to_process; ++k) {
                                            byte |= (bitset[pos + j + k] ? 1 : 0) << k;
                                        }
                                        simd_buffer[j/8] = byte;
                                    }
                                    bits = vld1q_u8(simd_buffer);
                                    vst1q_u8(simd_buffer, bits);
                                #endif

                                for (size_t j = 0; j < bits_to_process; ++j) {
                                    bins[bin][chunk_start + i + j] = (simd_buffer[j/8] >> (j%8)) & 1;
                                }
                            }
                        }
                    }
                }
            #else
                return partitionVector(bitset);
            #endif
            
            return bins;
        }

        boost::dynamic_bitset<> createBitVector(const std::vector<KeyType>& keys) {
            boost::dynamic_bitset<> bit_vector(range_);
            for (auto key : keys) {
                int offset = key - min_key_;
                bit_vector.set(offset); 
            }
            return bit_vector;
        }

        boost::dynamic_bitset<> createBitVectorSIMD(const std::vector<KeyType>& keys) {
            boost::dynamic_bitset<> bit_vector(range_);
            
            #ifdef __AVX2__
                __m256i min_key_vec = _mm256_set1_epi32(min_key_);
                for (size_t i = 0; i < keys.size(); i += 8) {
                    size_t remaining = keys.size() - i;
                    int temp_keys[8] = {0};

                    for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j) {
                        temp_keys[j] = keys[i + j];
                    }

                    __m256i key_vec = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(temp_keys));
                    __m256i offset_vec = _mm256_sub_epi32(key_vec, min_key_vec);

                    alignas(32) int offsets[8];
                    _mm256_store_si256(reinterpret_cast<__m256i*>(offsets), offset_vec);

                    for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j) {
                        if (offsets[j] >= 0 && offsets[j] < static_cast<int>(range_)) {
                            bit_vector.set(offsets[j]);
                        }
                    }
                }
            #elif defined(__ARM_NEON)
                int32x4_t min_key_vec = vdupq_n_s32(min_key_);
                for (size_t i = 0; i < keys.size(); i += 4) {
                    size_t remaining = keys.size() - i;
                    int temp_keys[4] = {0};

                    for (size_t j = 0; j < std::min<size_t>(4, remaining); ++j) {
                        temp_keys[j] = keys[i + j];
                    }

                    int32x4_t key_vec = vld1q_s32(temp_keys);
                    int32x4_t offset_vec = vsubq_s32(key_vec, min_key_vec);

                    int offsets[4];
                    vst1q_s32(offsets, offset_vec);

                    for (size_t j = 0; j < std::min<size_t>(4, remaining); ++j) {
                        if (offsets[j] >= 0 && offsets[j] < static_cast<int>(range_)) {
                            bit_vector.set(offsets[j]);
                        }
                    }
                }
            #else
                return createBitVector(keys);
            #endif

            return bit_vector;
        }


        static unsigned computeLog(uint32_t n, bool round = false) {
            assert(n);
            return 31 - __builtin_clz(n) + (round ? ((n & (n - 1)) != 0) : 0);
        }

        static unsigned computeLog(uint64_t n, bool round = false) {
            assert(n);
            return 63 - __builtin_clzl(n) + (round ? ((n & (n - 1)) != 0) : 0);
        }

        constexpr uint32_t setHighOrderBit(uint32_t value) {
            return value | (1 << 31);
        }

        constexpr bool isHighOrderBitSet(uint32_t value) {
            return value & (1 << 31);
        }

        constexpr uint32_t clearHighOrderBit(uint32_t value) {
            return value & ~(1 << 31);
        }

        static constexpr unsigned Terminal = 0xFFFFFFFF;

        KeyType min_key_;
        KeyType max_key_;
        const size_t num_bins_;
        const size_t log_num_bins_;
        const size_t max_error_;
        size_t range_;
        size_t next_inner_index = 0;
        size_t num_keys_;
        KeyType prev_key_;
        size_t shift_;

        std::vector<KeyType> keys_;
        std::vector<uint32_t> inner_nodes_;
        std::vector<uint32_t> leaf_nodes_;

};

