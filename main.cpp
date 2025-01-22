#include <immintrin.h> // Für AVX/AVX2 Intrinsics
#include <vector>
#include <iostream>
#include <bitset>
#include <cstdint>
#include <cmath>
#include <random>
#include <chrono>

#include "Builder.h"

char esc_char = 27;


int main() {
    std::vector<uint32_t> keys;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 250); // define the range

    size_t num_keys = 50;

    std::cout << "Generating " << num_keys << " random keys..." << std::endl;

    for(size_t n=0; n<num_keys; n++) {
        //check if the key is already in the vector
        uint32_t key = distr(gen);
        if (std::find(keys.begin(), keys.end(), key) != keys.end()) {
            n--;
            continue;
        }
        keys.push_back(key);
    }

    stable_sort(keys.begin(), keys.end());

    std::cout << esc_char << "[1m" << "The Keys: " << esc_char << "[0m";
    // print the sorted vector
    for (auto i : keys) {
        std::cout << i << " ";
    }   
    
    size_t maxerror = 4;
    size_t num_bins = 4;
    
    Builder<uint32_t> builder(keys, num_bins, maxerror);

    auto start = std::chrono::high_resolution_clock::now();
    HistTree<uint32_t> hist_tree = builder.build();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << std::endl << esc_char << "[1m" << "Time: " << esc_char << "[0m" << duration.count() << " microseconds" << std::endl;
    
    hist_tree.printVectors();

    hist_tree.visualize();
    
    std::cout << esc_char << "[1m" << "Input the key you want to search for or quit with q: " << esc_char << "[0m";
    std::string input;
    while (std::cin >> input) {
        if (input == "q") {
            break;
        }
        uint32_t key = std::stoi(input);
        SearchBound sb = hist_tree.getSearchBound(key);
        std::cout << esc_char << "[1m" << "Search Bound: " << esc_char << "[0m" << sb.start << " - " << sb.end << std::endl;
        std::cout << esc_char << "[1m" << "Input the key you want to search for or quit with q: " << esc_char << "[0m";
    }
    
    /*
    // same for insert/delete
    std::cout << esc_char << "[1m" << "Input the key you want to insert or quit with q: " << esc_char << "[0m";
    std::string input;
    while (std::cin >> input) {
        if (input == "q") {
            break;
        }
        uint32_t key = std::stoi(input);
        bool result = hist_tree.insert(key);
        std::cout << esc_char << "[1m" << "Insert Result: " << esc_char << "[0m" << result << std::endl;
        hist_tree.visualize();
        hist_tree.printVectors();
        std::cout << esc_char << "[1m" << "Input the key you want to insert or quit with q: " << esc_char << "[0m";
    }
    */
    
    return 0;
}
