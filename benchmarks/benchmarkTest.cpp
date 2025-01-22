#include <benchmark/benchmark.h>
#include <vector>

void ExampleFunction() {
    std::vector<int> vec(1000, 42);
    for (auto& v : vec) {
        v *= 2;
    }
}

static void BM_ExampleFunction(benchmark::State& state) {
    for (auto _ : state) {
        ExampleFunction();
    }
}

// register benchmark
BENCHMARK(BM_ExampleFunction);

BENCHMARK_MAIN();
