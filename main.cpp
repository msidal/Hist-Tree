#include <iostream>
#include <string>
#include <vector>
#include "Builder.h"
#include "HistTree.h"

void printUsage() {
    std::cout << "Commands:\n"
              << "insert <key>    - Insert a key\n"
              << "remove <key>    - Remove a key\n"
              << "search <key>    - Get search bounds for a key\n"
              << "quit            - Exit program\n";
}

int main() {
    // Tree configuration
    size_t num_bins, max_error, num_keys;
    std::vector<uint32_t> keys;
    
    std::cout << "Enter number of bins (power of 2): ";
    std::cin >> num_bins;
    
    std::cout << "Enter maximum error: ";
    std::cin >> max_error;
    
    std::cout << "Enter amount of keys: ";
    std::cin >> num_keys;
    
    std::cout << "Generating " << num_keys << " keys...\n";
    for (size_t i = 0; i < num_keys; i++) {
        keys.push_back(i);
    }

    
    // Build tree
    Builder<uint32_t> builder(keys, num_bins, max_error);
    auto tree = builder.build();
    
    // Command loop
    std::string command;
    uint32_t key;
    
    printUsage();
    
    while (true) {
        tree.visualize();
        std::cout << "\n> ";
        std::cin >> command;

        
        if (command == "quit") break;
        
        if (command == "insert") {
            std::cin >> key;
            if (tree.insert(key)) {
                std::cout << "Successfully inserted " << key << "\n";
            } else {
                std::cout << "Failed to insert " << key << "\n";
            }
        }
        else if (command == "remove") {
            std::cin >> key;
            if (tree.remove(key)) {
                std::cout << "Successfully removed " << key << "\n";
            } else {
                std::cout << "Failed to remove " << key << "\n";
            }
        }
        else if (command == "search") {
            std::cin >> key;
            auto bounds = tree.getSearchBound(key);
            std::cout << "Search Bound: [" << bounds.start << ", " << bounds.end << ")\n";
        }
        else {
            std::cout << "Unknown command\n";
            printUsage();
        }
    }
    
    return 0;
}