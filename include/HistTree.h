#pragma once

#include "omp.h"
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#ifdef USE_AVX2
#include <immintrin.h>
#endif

#include <common.h>

/**
 * @class HistTree
 * @brief A class representing a histogram tree for efficient range queries.
 *
 * @tparam KeyType The type of the keys stored in the tree.
 */
template <typename KeyType>
class HistTree
{
public:
    /**
     * @brief Default constructor.
     */
    HistTree() = default;

    /**
     * @brief Parameterized constructor.
     *
     * @param min_key Minimum key value.
     * @param max_key Maximum key value.
     * @param num_keys Number of keys.
     * @param num_bins Number of bins.
     * @param log_num_bins Logarithm of the number of bins.
     * @param max_error Maximum allowable error.
     * @param shift Shift value.
     * @param range Range of the keys.
     * @param inner_nodes Vector of inner nodes.
     * @param leaf_nodes Vector of leaf nodes.
     * @param bit_vector Bit vector representing the keys.
     */
    HistTree(KeyType min_key, KeyType max_key, size_t num_keys, size_t num_bins,
             size_t log_num_bins, size_t max_error, size_t shift, size_t range,
             std::vector<uint32_t> inner_nodes, std::vector<uint32_t> leaf_nodes,
             boost::dynamic_bitset<> bit_vector)
        : min_key_(min_key), max_key_(max_key), num_keys_(num_keys),
          num_bins_(num_bins), log_num_bins_(log_num_bins), max_error_(max_error),
          shift_(shift), range_(range), inner_nodes_(std::move(inner_nodes)),
          leaf_nodes_(std::move(leaf_nodes)),
          bit_vector_(std::move(bit_vector)) {};

    /**
     * @brief Get the search bound for a given key.
     *
     * @param key The key to search for.
     * @return SearchBound The search bound for the key.
     */
    SearchBound getSearchBound(const KeyType key) const
    {
        size_t begin = 0;
        try
        {
            begin = lookup(key);
        }
        catch (const std::out_of_range &e)
        {
            std::cerr << e.what() << std::endl;
            return SearchBound{0, 0};
        }

        // If the end of the range is greater than the number of keys
        const size_t end = (begin + max_error_ + 1 > num_keys_)
                               ? num_keys_
                               : (begin + max_error_ + 1);

        // `end` is exclusive
        return SearchBound{begin, end};
    }

    /**
     * @brief Lookup the first key equal to or greater than the given key.
     *
     * @param key The key to lookup.
     * @return size_t The position of the key.
     */
    size_t lookup(KeyType key) const
    {
        if (key < min_key_)
            return 0;
        if (key > max_key_)
            throw std::out_of_range("key out of range");

        size_t i, bin, pos = 0;
        key -= min_key_;

        // If root is a leaf node
        if (inner_nodes_.empty())
        {
            bin = key >> shift_;
            return std::accumulate(leaf_nodes_.begin(), leaf_nodes_.begin() + bin, 0);
        }

        uint32_t next, *node = const_cast<uint32_t *>(inner_nodes_.data());
        bool done = false;
        size_t width = shift_;

        do
        {
            // Find the bin
            bin = key >> width;

            // Accumulate the position
            for (i = 0; i < bin; i++)
            {
                pos += node[i];
            }

            if (done || node[num_bins_ + bin] == Terminal)
                return pos;

            done = isHighOrderBitSet(node[num_bins_ + bin]);
            next = clearHighOrderBit(node[num_bins_ + bin]);
            node = done ? const_cast<uint32_t *>(leaf_nodes_.data()) + next
                        : const_cast<uint32_t *>(inner_nodes_.data()) + next;

            // Update key and width
            key -= bin << width;
            width -= log_num_bins_;

        } while (1);
    }

    /**
     * @brief Insert a key into the tree.
     *
     * @param key The key to insert.
     * @return true If the key was successfully inserted.
     * @return false If the key already exists.
     */
    bool insert(KeyType key)
    {
        // Validate key and rebuild if necessary
        if (key < min_key_ || key > max_key_)
        {
            rebuild(key, RebuildContext::Insert);
            return true;
        }

        // Normalize key
        key -= min_key_;

        // Key already exists
        if (bit_vector_.at(key))
            return true; // idempotent

        // Insert key into bit vector
        bit_vector_.set(key);
        num_keys_++;

        size_t bin;

        // If root is a leaf node
        if (inner_nodes_.empty())
        {
            bin = key >> shift_;
            leaf_nodes_[bin]++;
            if (leaf_nodes_[bin] == max_error_)
                rebuild(key + min_key_, RebuildContext::Insert);
            return true;
        }

        uint32_t next, *node = inner_nodes_.data();
        bool done = false;
        size_t width = shift_;

        do
        {
            bin = key >> width;

            // Expand the tree
            if (++node[bin] == max_error_)
            {
                // could be optimized by "expanding in-place" but would use same algorithm as build
                rebuild(key + min_key_, RebuildContext::Insert);
                return true;
            }

            if (done || node[num_bins_ + bin] == Terminal)
                return true;

            done = isHighOrderBitSet(node[num_bins_ + bin]);
            next = clearHighOrderBit(node[num_bins_ + bin]);

            node = done ? leaf_nodes_.data() + next : inner_nodes_.data() + next;

            key -= bin << width;
            width -= log_num_bins_;

        } while (1);
    }

    /**
     * @brief Bulk insert multiple keys into the tree.
     *
     * @param keys The keys to insert.
     * @return std::vector<bool> A vector indicating the success of each insertion.
     */

    std::vector<bool> bulkInsert(const std::vector<KeyType> &keys)
    {
        std::vector<bool> success;
        for (const auto &key : keys)
        {
            success.push_back(insert(key));
        }
        return success;
    }

    /**
     * @brief Remove a key from the tree.
     *
     * @param key The key to remove.
     * @return true If the key was successfully removed.
     * @return false If the key does not exist.
     */
    bool remove(KeyType key)
    {
        if (key == min_key_ || key == max_key_)
        {
            rebuild(key, RebuildContext::Remove);
            return true;
        }

        // Validate key
        if (key < min_key_ || key > max_key_ || num_keys_ == 0 ||
            !bit_vector_.at(key - min_key_))
            return false;

        // Normalize key
        key -= min_key_;

        // Remove key from bit vector
        bit_vector_.reset(key);
        num_keys_--;

        size_t bin = 0;

        // If root is a leaf node
        if (inner_nodes_.empty())
        {
            bin = key >> shift_;
            leaf_nodes_[bin]--;
            return true;
        }

        // Traverse the tree
        uint32_t next, *node = inner_nodes_.data(), *parent = nullptr;
        size_t width = shift_, parent_bin = 0;
        bool is_root = true, is_leaf = false;
        do
        {
            bin = key >> width;
            node[bin]--;

            // If leaf node or terminal node reached
            if (is_leaf || node[num_bins_ + bin] == Terminal)
                return true;

            // Check if the node is an inner node or leaf node
            bool is_inner = std::any_of(node, node + num_bins_,
                                        [&](uint32_t n)
                                        {
                                            return n >= max_error_;
                                        });

            // If root becomes a leaf node
            if (is_root && !is_inner)
            {
                leaf_nodes_.assign(inner_nodes_.begin(),
                                   inner_nodes_.begin() + num_bins_);
                inner_nodes_.clear();
                break;
            }

            // If inner node becomes a leaf node
            if (!is_root && !is_inner)
            {
                // Create new leaf node
                uint32_t new_leaf_index = leaf_nodes_.size();
                for (size_t i = 0; i < num_bins_; i++)
                {
                    leaf_nodes_.push_back(node[i]);
                }

                // Update parent node
                parent[num_bins_ + parent_bin] = setHighOrderBit(new_leaf_index);

                // Cleanup the node and its children
                cleanupNode(node);

                break;
            }

            // If inner node becomes a terminal node
            if (node[bin] < max_error_ && is_inner)
            {
                if (isHighOrderBitSet(node[num_bins_ + bin]))
                {
                    uint32_t leaf_index = clearHighOrderBit(node[num_bins_ + bin]);
                    std::fill(leaf_nodes_.begin() + leaf_index,
                              leaf_nodes_.begin() + leaf_index + num_bins_, Filler);
                }
                else
                {
                    uint32_t *to_clean = &inner_nodes_[node[num_bins_ + bin]];
                    cleanupNode(to_clean);
                }

                // Update parent node
                node[num_bins_ + bin] = Terminal;
                break;
            }

            parent = node;
            parent_bin = bin;

            is_leaf = isHighOrderBitSet(node[num_bins_ + bin]);
            next = clearHighOrderBit(node[num_bins_ + bin]);
            node = is_leaf ? leaf_nodes_.data() + next : inner_nodes_.data() + next;

            key -= bin << width;
            width -= log_num_bins_;
            is_root = false;
        } while (1);

        return true;
    }

    /**
     * @brief Visualize the tree structure (recommended for small datasets).
     */
    void visualize() const
    {
        Visualize visualize(inner_nodes_, leaf_nodes_, num_bins_);
        visualize.exportToGraphviz("hist_tree.dot");
    };

    /**
     * @brief Get the size of the tree in bytes.
     *
     * @return size_t The size of the tree in bytes.
     */
    size_t getSize() const
    {
        // Size of the object itself
        size_t total_size = sizeof(*this);

        // Size of vector contents
        total_size += inner_nodes_.size() * sizeof(uint32_t);
        total_size += leaf_nodes_.size() * sizeof(uint32_t);

        // Vector capacity overhead
        total_size +=
            (inner_nodes_.capacity() - inner_nodes_.size()) * sizeof(uint32_t);
        total_size +=
            (leaf_nodes_.capacity() - leaf_nodes_.size()) * sizeof(uint32_t);

        // Size of vector objects themselves
        total_size += 3 * sizeof(void *); // for inner_nodes_
        total_size += 3 * sizeof(void *); // for leaf_nodes_

        // Size of boost::dynamic_bitset
        total_size += (bit_vector_.size() + 63) / 64 * sizeof(uint64_t);
        total_size += sizeof(size_t);

        return total_size;
    }

    /**
     * @brief To be implemented: Defragmantation.
     */
    void optimizeSpace()
    {
        inner_nodes_.shrink_to_fit();
        leaf_nodes_.shrink_to_fit();
    }

#ifdef TESTING
public:
#else
private:
#endif
    static constexpr unsigned Terminal = 0xFFFFFFFF; ///< Marks that the bin is terminal.
    static constexpr unsigned Filler = 0xFFFFFFFE;   ///< Marks freed space in the vectors.

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
    constexpr bool isHighOrderBitSet(uint32_t value) const
    {
        return value & (1 << 31);
    }

    /**
     * @brief Clear the high order bit of a value.
     *
     * @param value The value.
     * @return constexpr uint32_t The value with the high order bit cleared.
     */
    constexpr uint32_t clearHighOrderBit(uint32_t value) const
    {
        return value & ~(1 << 31);
    }

    /**
     * @brief Cleanup a node and its children.
     *
     * @param node The node to cleanup.
     */
    void cleanupNode(uint32_t *node)
    {
        // Queue of nodes to process
        std::queue<uint32_t> to_process;

        for (size_t i = num_bins_; i < num_bins_ * 2; i++)
        {
            if (node[i] != Terminal)
            {
                to_process.push(node[i]);
            }
        }

        std::fill(node, node + num_bins_ * 2, Filler);

        while (!to_process.empty())
        {
            uint32_t current = to_process.front();
            to_process.pop();

            if (isHighOrderBitSet(current))
            {
                // Leaf
                uint32_t leaf_index = clearHighOrderBit(current);
                std::fill(leaf_nodes_.begin() + leaf_index,
                          leaf_nodes_.begin() + leaf_index + num_bins_, Filler);
            }
            else
            {
                // Inner
                uint32_t *current_node = inner_nodes_.data() + current;
                for (size_t i = num_bins_; i < num_bins_ * 2; i++)
                {
                    if (current_node[i] != Terminal)
                    {
                        to_process.push(current_node[i]);
                    }
                }
                std::fill(current_node, current_node + num_bins_ * 2, Filler);
            }
        }
    }

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
     * @brief Setup the tree for rebuilding.
     *
     * @param key The key to rebuild around.
     * @param context The context of the rebuild (insert or remove).
     */
    void setupRebuild(KeyType key, RebuildContext context)
    {
        // regenerate the keys
        std::vector<KeyType> keys;
        keys.reserve(bit_vector_.count());
        for (size_t i = bit_vector_.find_first();
             i != boost::dynamic_bitset<>::npos; i = bit_vector_.find_next(i))
        {
            keys.push_back(i + min_key_);
        }

        if (context == RebuildContext::Insert)
        {
            keys.push_back(key);
            std::sort(keys.begin(), keys.end());
        }
        else
        {
            if (key == min_key_)
            {
                keys.erase(keys.begin());
            }
            else
            {
                keys.pop_back();
            }
        }

        // update the tree parameters
        min_key_ = keys.front();
        max_key_ = keys.back();

        auto log_range_ = computeLog(max_key_ - min_key_, true);
        shift_ = log_range_ - log_num_bins_;
        range_ = 1 << log_range_;
        bit_vector_ = createBitVectorSIMD(keys);
        num_keys_ = bit_vector_.count();
    }

    /**
     * @brief Partition the bitset into bins using SIMD.
     *
     * @param bitset The bitset to partition.
     * @return std::vector<boost::dynamic_bitset<>> Vector of bitsets representing the bins.
     * @details Check Builder.h for more information.
     */
    std::vector<boost::dynamic_bitset<>>
    partitionVectorSIMD(const boost::dynamic_bitset<> &bitset)
    {
        const size_t bin_size = bitset.size() / num_bins_;
        std::vector<boost::dynamic_bitset<>> bins(
            num_bins_, boost::dynamic_bitset<>(bin_size));

#ifdef USE_AVX2
        constexpr size_t simd_width = 256;
        constexpr size_t simd_bytes = simd_width / 8;
        const size_t L1_cache_size = 32768;
        const size_t chunk_size = (L1_cache_size / simd_bytes) * simd_width;

#pragma omp parallel
        {
            alignas(32) unsigned char simd_buffer[simd_bytes];

#pragma omp for schedule(dynamic)
            for (size_t bin = 0; bin < num_bins_; ++bin)
            {
                const size_t offset = bin * bin_size;

                for (size_t chunk_start = 0; chunk_start < bin_size;
                     chunk_start += chunk_size)
                {
                    const size_t current_chunk =
                        std::min(chunk_size, bin_size - chunk_start);

                    for (size_t i = 0; i < current_chunk; i += simd_width)
                    {
                        const size_t bits_to_process =
                            std::min(simd_width, current_chunk - i);
                        const size_t pos = offset + chunk_start + i;

                        __m256i bits = _mm256_setzero_si256();
                        for (size_t j = 0; j < bits_to_process; j += 8)
                        {
                            unsigned char byte = 0;
                            for (size_t k = 0; k < 8 && (j + k) < bits_to_process; ++k)
                            {
                                byte |= (bitset[pos + j + k] ? 1 : 0) << k;
                            }
                            simd_buffer[j / 8] = byte;
                        }
                        bits = _mm256_load_si256(
                            reinterpret_cast<const __m256i *>(simd_buffer));
                        _mm256_store_si256(reinterpret_cast<__m256i *>(simd_buffer), bits);

                        for (size_t j = 0; j < bits_to_process; ++j)
                        {
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
     * @brief Create a bit vector from the keys using SIMD.
     *
     * @param keys Vector of keys.
     * @return boost::dynamic_bitset<> The created bit vector.
     * @details Check Builder.h for more information.
     */
    boost::dynamic_bitset<>
    createBitVectorSIMD(const std::vector<KeyType> &keys)
    {
        boost::dynamic_bitset<> bit_vector(range_);

#ifdef USE_AVX2
        __m256i min_key_vec = _mm256_set1_epi32(min_key_);
        for (size_t i = 0; i < keys.size(); i += 8)
        {
            size_t remaining = keys.size() - i;
            int temp_keys[8] = {0};

            for (size_t j = 0; j < std::min<size_t>(8, remaining); ++j)
            {
                temp_keys[j] = keys[i + j];
            }

            __m256i key_vec =
                _mm256_loadu_si256(reinterpret_cast<const __m256i *>(temp_keys));
            __m256i offset_vec = _mm256_sub_epi32(key_vec, min_key_vec);

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
     * @brief Rebuild the tree.
     *
     * @param key The key to rebuild around.
     * @param context The context of the rebuild (insert or remove).
     * @details Same as Builder::build(). It is needed when the range changes. Could maybe be optimized to partially rebuild the tree.
     */
    void rebuild(KeyType key, RebuildContext context)
    {
        setupRebuild(key, context);

        size_t estimated_max_depth = static_cast<size_t>(
            std::ceil(std::log(static_cast<double>(num_keys_) / max_error_) /
                      std::log(num_bins_)));
        size_t estimated_inner_size =
            ((num_keys_ / max_error_) * (estimated_max_depth + 1) - 1)
            << (log_num_bins_ + 1);
        size_t estimated_leaf_size =
            static_cast<size_t>(std::pow(num_bins_, estimated_max_depth))
            << log_num_bins_;

        inner_nodes_.clear();
        leaf_nodes_.clear();
        inner_nodes_.resize(estimated_inner_size);
        leaf_nodes_.resize(estimated_leaf_size);

        auto bit_vector = bit_vector_;
        auto bins = partitionVectorSIMD(bit_vector);

        std::queue<std::tuple<std::vector<boost::dynamic_bitset<>>, size_t, size_t>>
            to_process;
        to_process.push({bins, 0, 0});

        size_t next_inner_index = 0;
        size_t next_leaf_index = 0;

        while (!to_process.empty())
        {
            auto [bins, current_index, bin_index] = to_process.front();
            to_process.pop();

            auto counts = countBinElements(bins);
            std::vector<uint32_t> child_count;

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
                for (size_t i = 0; i < num_bins_; ++i)
                {
                    inner_nodes_[current_index + i] = counts[i];
                }
                if (current_index == 0)
                    next_inner_index = num_bins_ << 1;

                for (size_t i = 0; i < num_bins_; ++i)
                {
                    if (counts[i] < max_error_)
                    {
                        inner_nodes_[current_index + num_bins_ + i] = Terminal;
                    }
                    else
                    {
                        auto child_bins = partitionVectorSIMD(bins[i]);
                        child_count = countBinElements(child_bins);
                        if (*std::max_element(child_count.begin(), child_count.end()) <
                            max_error_)
                        {
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
    }

    KeyType min_key_;           ///< Minimum key value.
    KeyType max_key_;           ///< Maximum key value.
    size_t num_keys_;           ///< Number of keys.
    const size_t num_bins_;     ///< Number of bins.
    const size_t log_num_bins_; ///< Logarithm of the number of bins.
    const size_t max_error_;    ///< Maximum allowable error.
    size_t shift_;              ///< Shift value.
    size_t range_;              ///< Range of the keys.

    std::vector<uint32_t> inner_nodes_; ///< Vector of inner nodes.
    std::vector<uint32_t> leaf_nodes_;  ///< Vector of leaf nodes.

    boost::dynamic_bitset<> bit_vector_; ///< Bit vector representing the keys (for updates).
};