#include <immintrin.h> // Für AVX/AVX2 Intrinsics
#include <vector>
#include <iostream>
#include <bitset>
#include <cstdint>
#include <cmath>
#include <random>

#include "Builder.h"


int main() {
     std::vector<uint32_t> keys;

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 10000); // define the range

    for(int n=0; n<4000; n++) {
        uint32_t key = distr(gen);
        while (std::find(keys.begin(), keys.end(), key) != keys.end())
        {
            key = distr(gen);
        }
        keys.push_back(key); // generate n numbers
        
    }

    stable_sort(keys.begin(), keys.end());

    size_t maxerror = 32;


    // erstelle builder
    Builder<uint32_t> builder(keys, 4, maxerror);

    HistTree<uint32_t> hist_tree = builder.build();

    builder.printVectors();


    return 0;
}
