#include <immintrin.h> // Für AVX/AVX2 Intrinsics
#include <vector>
#include <iostream>
#include <bitset>
#include <cstdint>
#include <cmath>
#include <random>
#include <chrono>

#include "Builder.h"


int main() {
     std::vector<uint32_t> keys;

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 110); // define the range

    for(int n=0; n<64; n++) {
        //check if the key is already in the vector
        uint32_t key = distr(gen);
        if (std::find(keys.begin(), keys.end(), key) != keys.end()) {
            n--;
            continue;
        }
        keys.push_back(key);   
    }

    stable_sort(keys.begin(), keys.end());

    // print the sorted vector
    for (auto i : keys) {
        std::cout << i << " ";
    }

    size_t maxerror = 2;

    
    Builder<uint32_t> builder(keys, 4, 2);

    auto start = std::chrono::high_resolution_clock::now();
    HistTree<uint32_t> hist_tree = builder.build();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Time: " << duration.count() << " microseconds" << std::endl;
    
    builder.printVectors();

    uint32_t cli_input = 0;
    std::cin >> cli_input;
    SearchBound sb = hist_tree.getSearchBound(cli_input);
    std::cout << "Search Bound: " << sb.start << " " << sb.end << std::endl;

    return 0;
}
