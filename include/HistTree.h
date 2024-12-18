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
            std::vector<uint32_t> inner_nodes, std::vector<uint32_t> leaf_nodes,
            std::vector<KeyType> keys)   :
            min_key_(min_key),
            max_key_(max_key),
            num_keys_(num_keys),
            num_bins_(num_bins),
            log_num_bins_(log_num_bins),
            max_error_(max_error),
            shift_(shift),
            inner_nodes_(std::move(inner_nodes)),
            leaf_nodes_(std::move(leaf_nodes)),
            keys_(std::move(keys)) {};

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
            for(i = 0; i < bin; i++) {
                pos += leaf_nodes_[i];
            }
            return pos;
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

    //TODO for the next three: I also need to update the keys, check if it is already in the tree, update min max
    bool insert(KeyType key) {
        return false;
    }


    //TODO is there a better way to do this?
    std::vector<bool> bulkInsert(const std::vector<KeyType>& keys) {
        std::vector<bool> success(keys.size(), true);
        for (const auto& key : keys) {
            success.push_back(insert(key));
        }
        return success;        
    }

    // TODO is it better to just build a new tree?
    bool remove(KeyType key) {
        // check if key is in range
        if (key < min_key_ || key > max_key_) return false;
        // check if tree is empty or key is not in the tree
        if (num_keys_ == 0 || keys_.count(key) == 0) return false;

        // remove key from keys
        keys_.erase(key);
        num_keys_--;

        size_t i, bin = 0;
        key -= min_key_;

        // if root is a leaf node
        if (inner_nodes_.empty()) {
            bin = key >> shift_;
            leaf_nodes_[bin]--;
            return true;
        }

        // traverse the tree
        uint32_t next, *node = inner_nodes_.data();
        size_t width = shift_;
        bool still_inner = false;

        /*
        do {
            bin = key >> width;
            node[bin]--;

            if (node[num_bins_ + bin] == Terminal) return true;

            for (i = 0; i < num_bins_; i++) {
                if(node[i] > max_error_) still_inner = true;
            }
        } while (1);
        */
        

        return true;
    }

    bool isEmpty() const {
        return num_keys_ == 0;
    }

    void visualize() const {
        Visualize visualize(inner_nodes_, leaf_nodes_, num_bins_);
        visualize.exportToGraphviz("hist_tree.dot");
    };

    // subject to deletion
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
 
    size_t getSize() const {
        return sizeof(*this) + inner_nodes_.size() * sizeof(uint32_t) + leaf_nodes_.size() * sizeof(uint32_t);
    }

    //TODO remove gaps in the vectors after a remove
    void optimize() {
        inner_nodes_.shrink_to_fit();
        leaf_nodes_.shrink_to_fit();
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

    // for inserts and removes
    std::vector<KeyType> keys_;

    KeyType min_key_;
    KeyType max_key_;
    size_t num_keys_;
    const size_t num_bins_;
    const size_t log_num_bins_;
    const size_t max_error_;
    const size_t shift_;

    // physical layout of the tree
    std::vector<uint32_t> inner_nodes_;
    std::vector<uint32_t> leaf_nodes_;

};
