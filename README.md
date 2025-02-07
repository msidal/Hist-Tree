# Hist-Tree

![License](https://img.shields.io/badge/license-MIT-green)
![C++](https://img.shields.io/badge/C++-17-blue.svg)

Hist-Tree [[1]](#1) is an efficient index competing with state-of-the-art Learned Indexes. It leverages histogram-based data partitioning with a tree structure for fast range queries and efficient memory usage.

## Features

- Support for 32-bit and 64-bit unsigned integers
- Range query lookup via SearchBound
- Dynamic updates through insert/remove operations  
- Visualization using Graphviz

## Requirements

Optional dependency for tree visualization:

```bash
sudo apt install graphviz
```

## Setup & Usage

Clone and run:

```bash
git clone https://github.com/msidal/hist-tree.git
cd hist-tree
chmod +x run.sh test.sh bench.sh
./run.sh
```

Available scripts:
- run.sh: Runs a simple CLI
- test.sh: Runs all unit tests
- bench.sh: Runs performance benchmarks

### Code Example

```cpp
#include "Builder.h"
#include "HistTree.h"

// Create tree with 4 bins and max error 2
std::vector<uint32_t> keys = {1, 2, 3, 4, 5};
Builder<uint32_t> builder(keys, 4, 2);
auto tree = builder.build();

// Search bounds
auto bounds = tree.getSearchBound(3);

// Dynamic updates
tree.insert(6);
tree.remove(1);
```

## References

<a id="1">[1]</a>
Andrew Crotty (2021).
Hist-Tree: Those Who Ignore It Are Doomed to Learn.
Conference on Innovative Data Systems Research. [Paper](https://api.semanticscholar.org/CorpusID:231400989)

## License

MIT License