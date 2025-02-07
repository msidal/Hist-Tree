#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <memory>
#include <iostream>
#include <numeric>
#include <boost/dynamic_bitset.hpp>

#ifdef __AVX2__
    #include <immintrin.h>
#elif defined(__ARM_NEON)
    #include <arm_neon.h>
#endif

#include <common.h>

template <typename KeyType>
class HistTree {
public:
    HistTree() = default; 

    HistTree(KeyType min_key, KeyType max_key, size_t num_keys, size_t num_bins, 
            size_t log_num_bins, size_t max_error, size_t shift, size_t range,
            std::vector<uint32_t> inner_nodes, std::vector<uint32_t> leaf_nodes,
            boost::dynamic_bitset<> bit_vector)   :
            min_key_(min_key),
            max_key_(max_key),
            num_keys_(num_keys),
            num_bins_(num_bins),
            log_num_bins_(log_num_bins),
            max_error_(max_error),
            shift_(shift),
            range_(range),
            inner_nodes_(std::move(inner_nodes)),
            leaf_nodes_(std::move(leaf_nodes)),
            bit_vector_(std::move(bit_vector)) {};

    SearchBound getSearchBound( const KeyType key) const {
        size_t begin = 0;
        try {
            begin = lookup(key);
        } catch (const std::out_of_range& e) {
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

    // lookup of first key equal or greater than key
    size_t lookup(KeyType key) const {
        if (key < min_key_) return 0;
        if (key > max_key_) throw std::out_of_range("key out of range");
        
        size_t i, bin, pos = 0;
        key -= min_key_;

        // if root is a leaf node
        if (inner_nodes_.empty()) {
            bin = key >> shift_;
            return std::accumulate(leaf_nodes_.begin(), leaf_nodes_.begin() + bin, 0);
        }

        uint32_t next, *node = const_cast<uint32_t*>(inner_nodes_.data());
        bool done = false;
        size_t width = shift_;

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

    bool insert(KeyType key) {
        // validate if key is in range
        if (key < min_key_ || key > max_key_) {
            rebuild(key, RebuildContext::Insert);
            return true;
        }
        
        // normalize key
        key -= min_key_;

        // validate if it already exists
        if (bit_vector_.at(key)) return true; // idempotent
        
        // insert key
        bit_vector_.set(key);
        num_keys_++;
        
        size_t bin;

        // if root is a leaf node
        if (inner_nodes_.empty()) {
            bin = key >> shift_;
            leaf_nodes_[bin]++;
            if (leaf_nodes_[bin] == max_error_) rebuild(key + min_key_, RebuildContext::Insert);
            return true;
        }

        uint32_t next, *node = inner_nodes_.data();
        bool done = false;
        size_t width = shift_;

        do {
            bin = key >> width;

            // expand the tree
            if (++node[bin] == max_error_) {
                // TODO hybrid approach
                rebuild(key + min_key_, RebuildContext::Insert);
                return true;
            }

            if (done || node[num_bins_ + bin] == Terminal) return true;

            done = isHighOrderBitSet(node[num_bins_ + bin]);
            next = clearHighOrderBit(node[num_bins_ + bin]);

            node = done ? leaf_nodes_.data() + next 
                        : inner_nodes_.data() + next;

            key -= bin << width;
            width -= log_num_bins_;

        } while (1);
    }

    std::vector<bool> bulkInsert(const std::vector<KeyType>& keys) {
        std::vector<bool> success;
        for (const auto& key : keys) {
            success.push_back(insert(key));
        }
        return success;        
    }

    bool remove(KeyType key) {
        if (key == min_key_ || key == max_key_) {
            rebuild(key, RebuildContext::Remove);
            return true;
        }

        // validate
        if (key < min_key_ || key > max_key_ || num_keys_ == 0 || !bit_vector_.at(key - min_key_)) return false;

        // normalize key
        key -= min_key_;

        // remove key from keys_ 
        bit_vector_.reset(key);
        num_keys_--;

        size_t bin = 0;

        // if root is a leaf node
        if (inner_nodes_.empty()) {
            bin = key >> shift_;
            leaf_nodes_[bin]--;
            return true; 
        }

        // traverse the tree
        uint32_t next, *node = inner_nodes_.data(), *parent = nullptr;
        size_t width = shift_, parent_bin = 0;
        bool is_root = true, is_leaf = false;
        do {
            // compute bin
            bin = key >> width;
            node[bin]--;

            // if leaf node is reached or the bin is terminal
            if (is_leaf || node[num_bins_ + bin] == Terminal) return true;
            
            // check if the node is an inner node after the decrement
            bool is_inner = std::any_of(node, node + num_bins_, [&](uint32_t n) {
                return n >= max_error_;
            });

            // if root becomes leaf node
            if (is_root && !is_inner) {
                leaf_nodes_.assign(inner_nodes_.begin(), inner_nodes_.begin() + num_bins_);
                inner_nodes_.clear();
                break;
            }

            if (!is_root && !is_inner) {
                // create new leaf node
                uint32_t new_leaf_index = leaf_nodes_.size();
                for(size_t i = 0; i < num_bins_; i++) {
                    leaf_nodes_.push_back(node[i]);
                }
        
                // update parent node
                parent[num_bins_ + parent_bin] = setHighOrderBit(new_leaf_index);

                // cleanup the node and its children
                cleanupNode(node);
               
                break;
            }

            if(node[bin] < max_error_ && is_inner) {
                if(isHighOrderBitSet(node[num_bins_ + bin])) {
                    uint32_t leaf_index = clearHighOrderBit(node[num_bins_ + bin]);
                    std::fill(leaf_nodes_.begin() + leaf_index, 
                            leaf_nodes_.begin() + leaf_index + num_bins_, 
                            Filler);
                } else {
                    uint32_t *to_clean = &inner_nodes_[node[num_bins_ + bin]];
                    cleanupNode(to_clean);
                }
                node[num_bins_ + bin] = Terminal;
                break;
            }

            parent = node;
            parent_bin = bin;

            is_leaf = isHighOrderBitSet(node[num_bins_ + bin]);
            next = clearHighOrderBit(node[num_bins_ + bin]);
            node = is_leaf ? leaf_nodes_.data() + next 
                        : inner_nodes_.data() + next;
            
            key -= bin << width;
            width -= log_num_bins_;
            is_root = false;
        } while (1);

        return true;
    }

    void printRange() const {
        std::cout << "Range: [" << min_key_ << ", " << max_key_ << "]" << std::endl;
    }

    // recommended to use only for debugging or small trees
    void visualize() const {
        Visualize visualize(inner_nodes_, leaf_nodes_, num_bins_);
        visualize.exportToGraphviz("hist_tree.dot");
    };

    // recommended to use only for debugging or small trees
    void printVectors() const {
        char esc_char = 27;
        std::cout << esc_char << "[1m"  << "Inner Nodes: " << esc_char << "[0m";
        for (const auto& node : inner_nodes_) {
            std::cout << node << " ";
        }
        std::cout << std::endl;

        std::cout << esc_char << "[1m"  << "Leaf Nodes: " << esc_char << "[0m";
        for (const auto& node : leaf_nodes_) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
 
    // Get size in bytes
    size_t getSize() const {
        // Size of the object itself
        size_t total_size = sizeof(*this);
        
        // Size of vector contents
        total_size += inner_nodes_.size() * sizeof(uint32_t);
        total_size += leaf_nodes_.size() * sizeof(uint32_t);
        
        // Vector capacity overhead
        total_size += (inner_nodes_.capacity() - inner_nodes_.size()) * sizeof(uint32_t);
        total_size += (leaf_nodes_.capacity() - leaf_nodes_.size()) * sizeof(uint32_t);
        
        // Size of vector objects themselves
        total_size += 3 * sizeof(void*); // for inner_nodes_
        total_size += 3 * sizeof(void*); // for leaf_nodes_
        
        // Size of boost::dynamic_bitset
        total_size += (bit_vector_.size() + 63) / 64 * sizeof(uint64_t);
        total_size += sizeof(size_t);
        
        return total_size;
    }

    //Future TODO: Implement a function to remove Filler elements from the vectors
    void optimizeSpace() {
        inner_nodes_.shrink_to_fit();
        leaf_nodes_.shrink_to_fit();
    }

#ifdef TESTING
    public:  
#else
    private: 
#endif
    static constexpr unsigned Terminal = 0xFFFFFFFF; // marks that the bin is terminal
    static constexpr unsigned Filler = 0xFFFFFFFE; //marks freed space in the vectors

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

    //helper to delete
    void cleanupNode(uint32_t* node) {
        std::queue<uint32_t> to_process;
        for (size_t i = num_bins_; i < num_bins_ * 2; i++) {
            if (node[i] != Terminal) {
                to_process.push(node[i]);
            }
        }

        std::fill(node, node + num_bins_ * 2, Filler);

        while(!to_process.empty()) {
            uint32_t current = to_process.front();
            to_process.pop();
            
            if (isHighOrderBitSet(current)) {
                uint32_t leaf_index = clearHighOrderBit(current);
                std::fill(leaf_nodes_.begin() + leaf_index, 
                        leaf_nodes_.begin() + leaf_index + num_bins_, 
                        Filler);
            } else {
                uint32_t* current_node = inner_nodes_.data() + current;
                for (size_t i = num_bins_; i < num_bins_ * 2; i++) {
                    if (current_node[i] != Terminal) {
                        to_process.push(current_node[i]);
                    }
                }
                std::fill(current_node, current_node + num_bins_ * 2, Filler);
            }
        }
    }

    // partition and count can later be updated/optimized just like in Builder.h

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


    void setupRebuild(KeyType key, RebuildContext context) {
        std::vector<KeyType> keys;
        keys.reserve(bit_vector_.count());
        for (size_t i = bit_vector_.find_first(); 
            i != boost::dynamic_bitset<>::npos;
            i = bit_vector_.find_next(i)) {
            keys.push_back(i + min_key_);
        }

        if (context == RebuildContext::Insert) {
            keys.push_back(key);
            std::sort(keys.begin(), keys.end());
        } else {
            if (key == min_key_) {
                keys.erase(keys.begin());
            } else  {
                keys.pop_back();
            } 
        }

        //update the tree parameters
        min_key_ = keys.front();
        max_key_ = keys.back();

        auto log_range_ = computeLog(max_key_ - min_key_, true);
        shift_ = log_range_ - log_num_bins_;
        range_ = 1 << log_range_;
        bit_vector_ = createBitVectorSIMD(keys);
        num_keys_ = bit_vector_.count();
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

    void rebuild(KeyType key, RebuildContext context) {
        setupRebuild(key, context);
        // estimated size of the tree
        size_t estimated_max_depth = static_cast<size_t>(std::ceil(std::log(static_cast<double>(num_keys_) / max_error_) / std::log(num_bins_)));
        size_t estimated_inner_size = ((num_keys_ / max_error_) * (estimated_max_depth + 1) - 1) << (log_num_bins_ + 1);
        size_t estimated_leaf_size = static_cast<size_t>(std::pow(num_bins_, estimated_max_depth)) << log_num_bins_;

        inner_nodes_.clear();
        leaf_nodes_.clear();
        inner_nodes_.resize(estimated_inner_size);
        leaf_nodes_.resize(estimated_leaf_size);

        auto bit_vector = bit_vector_;
        auto bins = partitionVectorSIMD(bit_vector); 
        
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
                        auto child_bins = partitionVectorSIMD(bins[i]);
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
    }

    KeyType min_key_;
    KeyType max_key_;
    size_t num_keys_;
    const size_t num_bins_;
    const size_t log_num_bins_;
    const size_t max_error_;
    size_t shift_;
    size_t range_;

    // physical layout of the tree
    std::vector<uint32_t> inner_nodes_;
    std::vector<uint32_t> leaf_nodes_;
     // for inserts and removes
    boost::dynamic_bitset<> bit_vector_;
};