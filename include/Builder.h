#pragma once

#include "HistTree.h"
// OpenMP for parallelism
#include "omp.h"
// faster than std::vector<bool> and more dynamic than std::bitset
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <queue>
#include <vector>

#ifdef USE_AVX2
#include <immintrin.h>
#endif

/**
 * @class Builder
 * @brief Builder class for constructing HistTree.
 *
 * @tparam KeyType Type of the keys.
 */
template <typename KeyType>
class Builder
{
public:
    /**
     * @brief Construct a new Builder object.
     *
     * @param keys Vector of keys.
     * @param num_bins Number of bins.
     * @param max_error Maximum allowable error.
     */
    Builder(std::vector<KeyType> keys, size_t num_bins, size_t max_error)
        : num_bins_(num_bins),
          log_num_bins_(computeLog(static_cast<uint64_t>(num_bins))),
          max_error_(max_error), num_keys_(keys.size()), keys_(std::move(keys))
    {
        assert(keys_.size() > 0 && "keys should not be empty");
        assert(num_bins_ <= num_keys_ &&
               "num_bins should not be greater than num_keys");
        assert(max_error_ <= num_keys_ &&
               "max_error should not be greater than num_keys");

        min_key_ = keys_.front();
        max_key_ = keys_.back();
        assert(min_key_ < max_key_);
        assert((num_bins_ & (num_bins_ - 1)) == 0 &&
               "num_bins should be a power of 2");
        assert(num_bins_ > 1 && "num_bins should be greater than 1");
        assert(max_error_ > 1 && "max_error should be greater than 1");

        auto log_range_ = computeLog(max_key_ - min_key_, true);
        assert(log_range_ >= log_num_bins_ &&
               "range should be greater than log_num_bins");
        shift_ = log_range_ - log_num_bins_;
        range_ = 1 << log_range_;
    }

    /**
     * @brief Build the HistTree.
     *
     * @return HistTree<KeyType> The constructed HistTree.
     */
    HistTree<KeyType> build()
    {
        // Estimate the size of the tree
        size_t estimated_max_depth = static_cast<size_t>(
            std::ceil(std::log(static_cast<double>(num_keys_) / max_error_) /
                      std::log(num_bins_)));
        size_t estimated_inner_size =
            ((num_keys_ / max_error_) * (estimated_max_depth + 1) - 1)
            << (log_num_bins_ + 1);
        size_t estimated_leaf_size =
            static_cast<size_t>(std::pow(num_bins_, estimated_max_depth))
            << log_num_bins_;

        inner_nodes_.resize(estimated_inner_size);
        leaf_nodes_.resize(estimated_leaf_size);

        // Create the bit vector and partition it into bins
        auto bit_vector = createBitVectorSIMD(keys_);
        auto bins = partitionVectorSIMD(bit_vector);

        // Queue of bins to process
        std::queue<std::tuple<std::vector<boost::dynamic_bitset<>>, size_t, size_t>>
            to_process;
        to_process.push({bins, 0, 0});

        size_t next_inner_index = 0;
        size_t next_leaf_index = 0;

        while (!to_process.empty())
        {
            auto [bins, current_index, bin_index] = to_process.front();
            to_process.pop();

            // Calculate the number of elements in each bin
            auto counts = countBinElements(bins);
            std::vector<uint32_t> child_count;

            // Check if the current node is a leaf
            if (*std::max_element(counts.begin(), counts.end()) < max_error_)
            {
                if (current_index > 0)
                {
                    inner_nodes_[current_index + num_bins_ + bin_index] =
                        setHighOrderBit(next_leaf_index);
                }
                for (size_t i = 0; i < num_bins_; ++i)
                {
                    leaf_nodes_[next_leaf_index + i] = counts[i];
                }
                next_leaf_index += num_bins_;
            }
            else
            {
                // Leaf
                for (size_t i = 0; i < num_bins_; ++i)
                {
                    inner_nodes_[current_index + i] = counts[i];
                }

                if (current_index == 0)
                    next_inner_index = num_bins_ << 1;

                // Process each bin
                for (size_t i = 0; i < num_bins_; ++i)
                {
                    // Check if the bin is terminal
                    if (counts[i] < max_error_)
                    {
                        inner_nodes_[current_index + num_bins_ + i] = Terminal;
                    }
                    else
                    {
                        // Partition the bin into child bins and count the elements
                        auto child_bins = partitionVectorSIMD(bins[i]);
                        child_count = countBinElements(child_bins);
                        if (*std::max_element(child_count.begin(), child_count.end()) <
                            max_error_)
                        {
                            // Create a leaf node
                            inner_nodes_[current_index + num_bins_ + i] =
                                setHighOrderBit(next_leaf_index);
                            for (size_t i = 0; i < counts.size(); ++i)
                            {
                                leaf_nodes_[next_leaf_index + i] = child_count[i];
                            }
                            next_leaf_index += num_bins_;
                        }
                        else
                        {
                            // Push the "inner node" child to the queue
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

        return HistTree<KeyType>(min_key_, max_key_, num_keys_, num_bins_,
                                 log_num_bins_, max_error_, shift_, range_,
                                 inner_nodes_, leaf_nodes_, bit_vector);
    }

#ifdef TESTING
public:
#else
private:
#endif
    /**
     * @brief Count the number of elements in each bin.
     *
     * @param bins Vector of bitsets representing the bins.
     * @return std::vector<uint32_t> Vector of counts for each bin.
     */
    std::vector<uint32_t>
    countBinElements(const std::vector<boost::dynamic_bitset<>> &bins)
    {
        std::vector<uint32_t> counts(num_bins_);
        for (size_t i = 0; i < num_bins_; ++i)
        {
            counts[i] = bins[i].count();
        }
        return counts;
    }

    /**
     * @brief Partition the bitset into bins.
     *
     * @param bitset The bitset to partition.
     * @return std::vector<boost::dynamic_bitset<>> Vector of bitsets representing the bins.
     */
    std::vector<boost::dynamic_bitset<>>
    partitionVector(const boost::dynamic_bitset<> &bitset)
    {
        size_t bin_size = bitset.size() / num_bins_;
        std::vector<boost::dynamic_bitset<>> result;
        result.reserve(num_bins_);
        size_t start_index = 0;

        for (size_t i = 0; i < num_bins_; ++i)
        {
            boost::dynamic_bitset<> bin(bin_size);
            for (size_t j = 0; j < bin_size; ++j)
            {
                bin[j] = bitset[start_index + j];
            }
            result.push_back(bin);
            start_index += bin_size;
        }
        return result;
    }

    /**
     * @brief Partition the bitset into bins using SIMD.
     *
     * @param bitset The bitset to partition.
     * @return std::vector<boost::dynamic_bitset<>> Vector of bitsets representing the bins.
     * @details This optimized version of partitionVector() uses SIMD instructions (AVX2)
     * to reduce its overhead during the build process (see docs/FlameGraphBuild.svg),
     * making it particularly efficient for large key sets.
     */
    std::vector<boost::dynamic_bitset<>>
    partitionVectorSIMD(const boost::dynamic_bitset<> &bitset)
    {
        const size_t bin_size = bitset.size() / num_bins_;
        std::vector<boost::dynamic_bitset<>> bins(
            num_bins_, boost::dynamic_bitset<>(bin_size));

#ifdef USE_AVX2
        // SIMD parameters
        constexpr size_t simd_width = 256;
        constexpr size_t simd_bytes = simd_width / 8;
        const size_t L1_cache_size = 32768;
        const size_t chunk_size = (L1_cache_size / simd_bytes) * simd_width;

#pragma omp parallel
        {
            // SIMD buffer
            alignas(32) unsigned char simd_buffer[simd_bytes];

#pragma omp for schedule(dynamic)
            for (size_t bin = 0; bin < num_bins_; ++bin)
            {
                const size_t offset = bin * bin_size;

                // Process the bitset in chunks
                for (size_t chunk_start = 0; chunk_start < bin_size;
                     chunk_start += chunk_size)
                {
                    const size_t current_chunk =
                        std::min(chunk_size, bin_size - chunk_start);

                    for (size_t i = 0; i < current_chunk; i += simd_width)
                    {
                        const size_t bits_to_process =
                            std::min(simd_width, current_chunk - i);
                        // Absolute position in the bitset
                        const size_t pos = offset + chunk_start + i;

                        // Initialize the SIMD buffer
                        __m256i bits = _mm256_setzero_si256();

                        // Put respectively 8 bits into the SIMD buffer
                        for (size_t j = 0; j < bits_to_process; j += 8)
                        {
                            unsigned char byte = 0;
                            // Construct the byte
                            for (size_t k = 0; k < 8 && (j + k) < bits_to_process; ++k)
                            {
                                // Set the k-th bit of the byte if the bitset bit is set
                                byte |= (bitset[pos + j + k] ? 1 : 0) << k;
                            }
                            simd_buffer[j / 8] = byte;
                        }
                        // Load the SIMD buffer into the SIMD register
                        bits = _mm256_load_si256(
                            reinterpret_cast<const __m256i *>(simd_buffer));
                        // Store the SIMD register into the SIMD buffer
                        _mm256_store_si256(reinterpret_cast<__m256i *>(simd_buffer), bits);

                        // Unpack the SIMD buffer into the bin
                        for (size_t j = 0; j < bits_to_process; ++j)
                        {
                            // Extract the j-th bit from the SIMD buffer
                            bins[bin][chunk_start + i + j] =
                                (simd_buffer[j / 8] >> (j % 8)) & 1;
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

    /**
     * @brief Create a bit vector from the keys.
     *
     * @param keys Vector of keys.
     * @return boost::dynamic_bitset<> The created bit vector.
     */
    boost::dynamic_bitset<> createBitVector(const std::vector<KeyType> &keys)
    {
        boost::dynamic_bitset<> bit_vector(range_);
        for (auto key : keys)
        {
            int offset = key - min_key_;
            bit_vector.set(offset);
        }
        return bit_vector;
    }

    /**
     * @brief Create a bit vector from the keys using SIMD.
     *
     * @param keys Vector of keys.
     * @return boost::dynamic_bitset<> The created bit vector.
     * @details This SIMD-optimized version processes 8 keys in parallel using AVX2 instructions,
     * making it particularly efficient for large key sets.
     */
    boost::dynamic_bitset<>
    createBitVectorSIMD(const std::vector<KeyType> &keys)
    {
        boost::dynamic_bitset<> bit_vector(range_);

#ifdef USE_AVX2
        // Broadcast the minimum key
        __m256i min_key_vec = _mm256_set1_epi32(min_key_);

        // Process the keys in chunks of 8
        for (size_t i = 0; i < keys.size(); i += 8)
        {
            size_t remaining = keys.size() - i;
            int temp_keys[8] = {0};

            for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j)
            {
                temp_keys[j] = keys[i + j];
            }

            // Load the keys into a SIMD register
            __m256i key_vec =
                _mm256_loadu_si256(reinterpret_cast<const __m256i *>(temp_keys));
            // Compute the offset of the keys
            __m256i offset_vec = _mm256_sub_epi32(key_vec, min_key_vec);

            // Store the offsets into an array
            alignas(32) int offsets[8];
            _mm256_store_si256(reinterpret_cast<__m256i *>(offsets), offset_vec);

            for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j)
            {
                if (offsets[j] >= 0 && offsets[j] < static_cast<int>(range_))
                {
                    bit_vector.set(offsets[j]);
                }
            }
        }
#else
        return createBitVector(keys);
#endif

        return bit_vector;
    }

    /**
     * @brief Compute the logarithm base 2 of a number.
     *
     * @param n The number.
     * @param round Whether to round up.
     * @return unsigned The computed logarithm.
     */
    static unsigned computeLog(uint32_t n, bool round = false)
    {
        assert(n);
        return 31 - __builtin_clz(n) + (round ? ((n & (n - 1)) != 0) : 0);
    }

    /**
     * @brief Compute the logarithm base 2 of a number.
     *
     * @param n The number.
     * @param round Whether to round up.
     * @return unsigned The computed logarithm.
     */
    static unsigned computeLog(uint64_t n, bool round = false)
    {
        assert(n);
        return 63 - __builtin_clzl(n) + (round ? ((n & (n - 1)) != 0) : 0);
    }

    /**
     * @brief Set the high order bit of a value.
     *
     * @param value The value.
     * @return constexpr uint32_t The value with the high order bit set.
     */
    constexpr uint32_t setHighOrderBit(uint32_t value)
    {
        return value | (1 << 31);
    }

    /**
     * @brief Check if the high order bit is set.
     *
     * @param value The value.
     * @return constexpr bool True if the high order bit is set, false otherwise.
     */
    constexpr bool isHighOrderBitSet(uint32_t value)
    {
        return value & (1 << 31);
    }

    /**
     * @brief Clear the high order bit of a value.
     *
     * @param value The value.
     * @return constexpr uint32_t The value with the high order bit cleared.
     */
    constexpr uint32_t clearHighOrderBit(uint32_t value)
    {
        return value & ~(1 << 31);
    }

    static constexpr unsigned Terminal = 0xFFFFFFFF; ///< Terminal value for inner nodes.

    KeyType min_key_;           ///< Minimum key.
    KeyType max_key_;           ///< Maximum key.
    const size_t num_bins_;     ///< Number of bins.
    const size_t log_num_bins_; ///< Logarithm of the number of bins.
    const size_t max_error_;    ///< Maximum allowable error.
    size_t range_;              ///< Range of the keys.
    size_t num_keys_;           ///< Number of keys.
    size_t shift_;              ///< Shift value.

    std::vector<KeyType> keys_; ///< Vector of keys.

    // physical representation of the tree
    std::vector<uint32_t> inner_nodes_; ///< Inner nodes.
    std::vector<uint32_t> leaf_nodes_;  ///< Leaf nodes.
};
