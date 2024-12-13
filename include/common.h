#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <queue>

struct SearchBound
{
    size_t start;
    size_t end; // exclusive
};


class Visualize {
public:
    Visualize(const std::vector<uint32_t>& inner_nodes, 
              const std::vector<uint32_t>& leaf_nodes, 
              size_t num_bins,
              size_t max_error) 
        : inner_nodes_(inner_nodes), 
          leaf_nodes_(leaf_nodes), 
          num_bins_(num_bins),
          max_error_(max_error) {}

    // creates a folder "graphs" and a file "filename" in it
    static void createGraphsFolderAndFile(const std::string& filename) {
        auto rootPath = std::filesystem::current_path().parent_path(); 
        std::filesystem::path folderPath = rootPath / "graphs";
        std::filesystem::path filePath = folderPath / filename;

        if (!std::filesystem::exists(folderPath)) {
            if (std::filesystem::create_directory(folderPath)) {
                std::cout << "Folder 'graphs' was created.\n";
            } else {
                std::cerr << "Error creating folder 'graphs'.\n";
                return;
            }
        }

        if (!std::filesystem::exists(filePath)) {
            std::ofstream file(filePath);
            if (file) {
                std::cout << "File '" << filename << "' was created in 'graphs'.\n";
            } else {
                std::cerr << "Error creating file '" << filename << "'.\n";
            }
        } else {
            std::cout << "File '" << filename << "' already exists in 'graphs'.\n";
        }
    }


    void exportToGraphviz(const std::string& filename) const {
        createGraphsFolderAndFile(filename);
        auto rootPath = std::filesystem::current_path().parent_path(); 
        std::filesystem::path folderPath = rootPath / "graphs";
        std::filesystem::path filePath = folderPath / filename;

        std::ofstream out(filePath);
        if (!out) {
            std::cerr << "Error opening the file: " << filePath << std::endl;
            return;
        }

        // configure the graph
        out << "digraph HistTree {\n";
        out << "  node [shape=record, fontname=Helvetica];\n";  // record shape for nodes
        out << "  rankdir=TB;\n"; // top to bottom
        out << "  ranksep=1.0;\n";
        out << "  nodesep=0.5;\n";
        out << "  size=\"12,24\";\n";

        // start with the root node
        if (!inner_nodes_.empty()) {
            exportNode(out, 0, true);
        }

        out << "}\n";
        out.close();
        std::cout << "Graphviz file has been created." << filePath << std::endl;
        std::string command = "dot -Tpng " + filePath.string() + " -o " + folderPath.string() + "/" + replaceDotWithPng(filename);
        system(command.c_str());
        filePath.replace_extension(".png");
        std::cout << "Graphviz file has been converted to PNG in folder 'graphs'." << filePath << std::endl;
    }


private:
    // helper function to export a node
    void exportNode(std::ofstream& out, size_t index, bool isInner) const {
        const std::vector<uint32_t>& nodes = isInner ? inner_nodes_ : leaf_nodes_;

        // create record label for the node
        std::string label = "";
        for (size_t i = 0; i < num_bins_; ++i) {
            label += std::to_string(nodes[index + i]) + (i < num_bins_ - 1 ? "|" : "");
        }
        
        out << "  node" << index << " [label=\"" << label << "}\", shape=record];\n";

        // export all children
        for (size_t i = 0; i < num_bins_; ++i) {
            size_t child_index = 0;
            child_index = index + num_bins_+ i;

            if (nodes[child_index] == Terminal) {
                continue; // ignore terminal nodes
            }

            // leaf node
            if (isHighOrderBitSet(nodes[child_index])) {
                size_t leaf_index = clearHighOrderBit(nodes[child_index]);
                out << "  node" << index << " -> leaf_" << leaf_index << " [label=\"" << i << "\"];\n";
                exportLeafNode(out, leaf_index);  
            // inner node
            } else {
                size_t next_index = nodes[child_index];
                out << "  node" << index << " -> node" << next_index << " [label=\"" << i << "\"];\n";
                exportNode(out, next_index, true);
            }
        }
    }

    // helper function to export a leaf node
    void exportLeafNode(std::ofstream& out, size_t leaf_index) const {
        const uint32_t value = leaf_nodes_[leaf_index];
        std::string label = std::to_string(value);
        
        std::string fields = "";
        for (size_t i = 0; i < num_bins_; ++i) {
            fields += std::to_string(leaf_nodes_[leaf_index + i]) + (i < num_bins_ - 1 ? "|" : "");
        }

        out << "  leaf_" << leaf_index << " [label=\"{" << fields << "}\", shape=record];\n";
    }

    constexpr static uint32_t Terminal = 0xFFFFFFFF;
    

    // helper functions for high order bit 
    static constexpr bool isHighOrderBitSet(uint32_t value) {
        return value & (1u << 31);
    }

    static constexpr uint32_t clearHighOrderBit(uint32_t value) {
        return value & ~(1u << 31);
    }

    // helper function to replace .dot with .png
    std::string replaceDotWithPng(const std::string& filename) const {
        std::string newFilename = filename;
        size_t pos = newFilename.rfind(".dot"); // Suche nach der Endung .dot
        if (pos != std::string::npos) {
            newFilename.replace(pos, 4, ".png"); // Ersetze .dot durch .png
        }
        return newFilename;
    }


    const std::vector<uint32_t>& inner_nodes_;
    const std::vector<uint32_t>& leaf_nodes_;
    size_t num_bins_;
    size_t max_error_;
};

