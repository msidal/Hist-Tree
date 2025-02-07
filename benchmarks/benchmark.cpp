#define TESTING
#include <benchmark/benchmark.h>
#include <random>
#include <vector>
#include <algorithm>
#include <limits>
#include <memory>
#include "omp.h"
#include "Builder.h"
#include "HistTree.h"
#include <set>

/*
#include <set> // für std::set, welches typischerweise als Rot-Schwarz-Baum implementiert ist

// BTree Benchmark für Vergleich
static void BM_BTreeBuildPerformance(benchmark::State& state) {
    const size_t dataSize = state.range(0);
    
    auto data = generateSortedData<uint32_t>(dataSize);
    double totalTreeSize = 0.0;

    for (auto _ : state) {
        state.PauseTiming();
        std::set<uint32_t> btree; // oder alternative BTree-Implementierung
        state.ResumeTiming();

        for (const auto& key : data) {
            btree.insert(key);
        }

        state.PauseTiming();
        totalTreeSize = static_cast<double>(btree.size() * sizeof(uint32_t) * 3) / (1024.0 * 1024.0); // Schätzung für BTree Overhead
        state.ResumeTiming();
    }
   
    state.counters["TotalTreeSize_MB"] = totalTreeSize;
    state.SetComplexityN(dataSize);
    state.SetLabel("BTree - Size: " + std::to_string(dataSize));
}

// Dann die Registrierung:
BENCHMARK(BM_BTreeBuildPerformance)
    ->Args({2000000})
    ->Args({20000000})
    ->Args({200000000})
    ->Unit(benchmark::kMicrosecond)
    ->Complexity(benchmark::oNLogN)
    ->UseRealTime();
*/

template <typename KeyType>
std::vector<KeyType> generateSortedData(size_t size) {
    std::vector<KeyType> data(size);
    #pragma omp parallel for
    for (size_t i = 0; i < size; i++) {
        // Clustered around center with controlled spread and guaranteed uniqueness
        data[i] = static_cast<KeyType>(size/2 + i + std::sin(i) * size * 0.05);
    }
    std::sort(data.begin(), data.end());
    return data;
}

std::string formatLabel(size_t size, size_t bins, size_t error) {
   return "Size: " + std::to_string(size) + 
          ", Bins: " + std::to_string(bins) + 
          ", MaxError: " + std::to_string(error);
}

static void BM_HistTreeBuildPerformance(benchmark::State& state) {
    const size_t dataSize = state.range(0);
    const size_t numBins = state.range(1);
    const size_t maxError = state.range(2);
    
    auto data = generateSortedData<uint32_t>(dataSize);
    double totalTreeSize = 0.0;

    for (auto _ : state) {
        try {
            state.PauseTiming();
            Builder<uint32_t> builder(data, numBins, maxError);
            state.ResumeTiming();

            auto tree = builder.build();

            state.PauseTiming();
            totalTreeSize = static_cast<double>(tree.getSize()) / (1024.0 * 1024.0);
            state.ResumeTiming();
       } catch (const std::exception& e) {
            state.SkipWithError(("Build error: " + std::string(e.what())).c_str());
            break;
        }
    }
   
    // Use MB in counters
    state.counters["TotalTreeSize_MB"] = totalTreeSize;
    state.SetComplexityN(dataSize);
    state.SetLabel(formatLabel(dataSize, numBins, maxError));
}

static void BM_HistTreeGetSearchBound(benchmark::State& state) {
    const size_t dataSize = state.range(0);
    const size_t numBins = state.range(1);
    const size_t maxError = state.range(2);
    
    auto data = generateSortedData<uint32_t>(dataSize);
    Builder<uint32_t> builder(data, numBins, maxError);
    auto tree = builder.build();
    
    for (auto _ : state) {
        state.PauseTiming();
        uint32_t key = data[rand() % data.size()];
        state.ResumeTiming();
        benchmark::DoNotOptimize(tree.getSearchBound(key));
    }
    
    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(dataSize);
    state.SetLabel(formatLabel(dataSize, numBins, maxError));
}

static void BM_HistTreeRemove(benchmark::State& state) {
    const size_t dataSize = state.range(0);
    const size_t numBins = state.range(1);
    const size_t maxError = state.range(2);
    
    auto data = generateSortedData<uint32_t>(dataSize);
    uint32_t key = data[dataSize/2];
    
    auto tree = Builder<uint32_t>(data, numBins, maxError).build();
    
    for (auto _ : state) {
        benchmark::DoNotOptimize(tree.remove(key));
    }
    
    state.SetItemsProcessed(1);
    state.SetComplexityN(dataSize);
    state.SetLabel(formatLabel(dataSize, numBins, maxError));
}

static void BM_HistTreeInsert(benchmark::State& state) {
    const size_t dataSize = state.range(0);
    const size_t numBins = state.range(1);
    const size_t maxError = state.range(2);
    
    auto data = generateSortedData<uint32_t>(dataSize);
    uint32_t key = data[dataSize/2];
    
    auto tree = Builder<uint32_t>(data, numBins, maxError).build();
    
    for (auto _ : state) {
        benchmark::DoNotOptimize(tree.insert(key + 1));
    }
    
    state.SetItemsProcessed(1);
    state.SetComplexityN(dataSize);
    state.SetLabel(formatLabel(dataSize, numBins, maxError));
}

/*
BENCHMARK(BM_HistTreeBuildPerformance)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Unit(benchmark::kMicrosecond)
    ->Complexity(benchmark::oNLogN)
    ->UseRealTime();

BENCHMARK(BM_HistTreeInsert)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Iterations(1)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(BM_HistTreeGetSearchBound)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(BM_HistTreeRemove)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Iterations(1)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
*/

/*
template<typename KeyType>
struct BinarySearchIndex {
    std::vector<KeyType> data;
    
    BinarySearchIndex(const std::vector<KeyType>& input) : data(input) {}
    
    size_t binary_search(KeyType target) const {
        size_t left = 0;
        size_t right = data.size() - 1;
        
        while (left <= right) {
            size_t mid = left + (right - left) / 2;
            
            if (data[mid] == target) {
                return mid;
            }
            
            if (data[mid] < target) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        return left; // Return insert position if not found
    }

    // Analog zu getSearchBound im HistTree
    std::pair<size_t, size_t> getSearchBound(KeyType key, size_t max_error) const {
        size_t pos = binary_search(key);
        size_t end = std::min(pos + max_error + 1, data.size());
        return {pos, end};
    }

    size_t get_memory_usage() const {
        return sizeof(BinarySearchIndex) + data.size() * sizeof(KeyType);
    }
};


static void BM_BinarySearchGetBound(benchmark::State& state) {
    const size_t dataSize = state.range(0);
    const size_t maxError = state.range(2);
    
    auto data = generateSortedData<uint32_t>(dataSize);
    BinarySearchIndex<uint32_t> index(data);
    
    for (auto _ : state) {
        state.PauseTiming();
        uint32_t key = data[rand() % data.size()];
        state.ResumeTiming();
        
        auto bound = index.getSearchBound(key, maxError);
        benchmark::DoNotOptimize(bound);
    }
    
    state.SetItemsProcessed(state.iterations());
    state.SetComplexityN(dataSize);
    state.SetLabel("Binary Search - Size: " + std::to_string(dataSize));
}


BENCHMARK(BM_BinarySearchGetBound)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {64},  // Dummy für Vergleichbarkeit
        {2048} // maxError
    })
    ->Unit(benchmark::kNanosecond)
    ->Complexity();
*/

BENCHMARK(BM_HistTreeBuildPerformance)
    ->ArgsProduct({
        {200000000},
        {64},
        {2048}
    })
    ->Unit(benchmark::kMicrosecond)
    ->Complexity(benchmark::oNLogN)
    ->UseRealTime();

BENCHMARK(BM_HistTreeInsert)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Iterations(1)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(BM_HistTreeGetSearchBound)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Unit(benchmark::kNanosecond)
    ->Complexity();

BENCHMARK(BM_HistTreeRemove)
    ->ArgsProduct({
        {2000000, 20000000, 200000000},
        {32, 64, 128},
        {1024, 2048, 8192}
    })
    ->Iterations(1)
    ->Unit(benchmark::kNanosecond)
    ->Complexity();


BENCHMARK_MAIN();