#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

/**
 * @struct SearchBound
 * @brief Represents the search boundaries within a histogram tree.
 */
struct SearchBound
{
    size_t start; ///< The starting index of the search boundary.
    size_t end;   ///< The ending index of the search boundary (exclusive).
};

/**
 * @enum RebuildContext
 * @brief Enumerates the contexts in which a histogram tree might be rebuilt.
 */
enum RebuildContext
{
    Insert, ///< Context for inserting nodes.
    Remove  ///< Context for removing nodes.
};

/**
 * @class Visualize
 * @brief A class to visualize histogram trees and export them to Graphviz format.
 */
class Visualize
{
public:
    /**
     * @brief Constructor for the Visualize class.
     * @param inner_nodes A vector of inner node values.
     * @param leaf_nodes A vector of leaf node values.
     * @param num_bins The number of bins for each node.
     */
    Visualize(const std::vector<uint32_t> &inner_nodes,
              const std::vector<uint32_t> &leaf_nodes, size_t num_bins)
        : inner_nodes_(inner_nodes), leaf_nodes_(leaf_nodes),
          num_bins_(num_bins) {}

    /**
     * @brief Creates a folder "graphs" and a file with the specified filename in it.
     * @param filename The name of the file to be created.
     */
    static void createGraphsFolderAndFile(const std::string &filename)
    {
        auto rootPath = std::filesystem::current_path().parent_path();
        std::filesystem::path folderPath = rootPath / "graphs";
        std::filesystem::path filePath = folderPath / filename;

        if (!std::filesystem::exists(folderPath))
        {
            if (std::filesystem::create_directory(folderPath))
            {
                std::cout << "Folder 'graphs' was created.\n";
            }
            else
            {
                std::cerr << "Error creating folder 'graphs'.\n";
                return;
            }
        }

        if (!std::filesystem::exists(filePath))
        {
            std::ofstream file(filePath);
            if (file)
            {
                std::cout << "File '" << filename << "' was created in 'graphs'.\n";
            }
            else
            {
                std::cerr << "Error creating file '" << filename << "'.\n";
            }
        }
        else
        {
            std::cout << "File '" << filename << "' already exists in 'graphs'.\n";
        }
    }

    /**
     * @brief Exports the histogram tree data to a Graphviz file and converts it to a PNG image.
     * @param filename The name of the Graphviz file to be created.
     */
    void exportToGraphviz(const std::string &filename) const
    {
        createGraphsFolderAndFile(filename);
        auto rootPath = std::filesystem::current_path().parent_path();
        std::filesystem::path folderPath = rootPath / "graphs";
        std::filesystem::path filePath = folderPath / filename;

        std::ofstream out(filePath);
        if (!out)
        {
            std::cerr << "Error opening the file: " << filePath << std::endl;
            return;
        }

        // Configure the graph
        out << "digraph HistTree {\n";
        out << "  node [shape=record, fontname=Helvetica];\n"; // Record shape for nodes
        out << "  rankdir=TB;\n";                              // Top to bottom
        out << "  ranksep=1.0;\n";
        out << "  nodesep=0.5;\n";

        // Start with the root node

        // If the tree is empty, create a root node with no children
        if (inner_nodes_.empty())
        {
            std::string label = "";
            for (size_t i = 0; i < num_bins_; ++i)
            {
                label +=
                    std::to_string(leaf_nodes_[i]) + (i < num_bins_ - 1 ? "|" : "");
            }
            out << "  root [label=\"" << label << "\", shape=record];\n";
        }
        else
        {
            exportNode(out, 0, true);
        }

        out << "}\n";
        out.close();
        std::cout << "Graphviz file has been created." << filePath << std::endl;
        std::string command = "dot -Tpng " + filePath.string() + " -o " +
                              folderPath.string() + "/" +
                              replaceDotWithPng(filename);
        auto success = system(command.c_str());
        (void)success;
        filePath.replace_extension(".png");
        std::cout << "Graphviz file has been converted to PNG in folder 'graphs'."
                  << filePath << std::endl;
    }

private:
    /**
     * @brief Helper function to export a node to the Graphviz file.
     * @param out The output file stream.
     * @param index The index of the node to be exported.
     * @param isInner A flag indicating whether the node is an inner node.
     */
    void exportNode(std::ofstream &out, size_t index, bool isInner) const
    {
        const std::vector<uint32_t> &nodes = isInner ? inner_nodes_ : leaf_nodes_;

        // Create record label for the node
        std::string label = "";
        for (size_t i = 0; i < num_bins_; ++i)
        {
            label +=
                std::to_string(nodes[index + i]) + (i < num_bins_ - 1 ? "|" : "");
        }

        out << "  node" << index << " [label=\"" << label
            << "}\", shape=record];\n";

        // Export all children
        for (size_t i = 0; i < num_bins_; ++i)
        {
            size_t child_index = 0;
            child_index = index + num_bins_ + i;

            if (nodes[child_index] == Terminal || nodes[child_index] == Filler)
            {
                continue; // ignore terminal nodes
            }

            if (isHighOrderBitSet(nodes[child_index]))
            {
                // Leaf node
                size_t leaf_index = clearHighOrderBit(nodes[child_index]);
                out << "  node" << index << " -> leaf_" << leaf_index << " [label=\""
                    << i << "\"];\n";
                exportLeafNode(out, leaf_index);
                
            }
            else
            {
                // Inner node
                size_t next_index = nodes[child_index];
                out << "  node" << index << " -> node" << next_index << " [label=\""
                    << i << "\"];\n";
                exportNode(out, next_index, true);
            }
        }
    }

    /**
     * @brief Helper function to export a leaf node to the Graphviz file.
     * @param out The output file stream.
     * @param leaf_index The index of the leaf node to be exported.
     */
    void exportLeafNode(std::ofstream &out, size_t leaf_index) const
    {
        const uint32_t value = leaf_nodes_[leaf_index];
        std::string label = std::to_string(value);

        std::string fields = "";
        for (size_t i = 0; i < num_bins_; ++i)
        {
            fields += std::to_string(leaf_nodes_[leaf_index + i]) +
                      (i < num_bins_ - 1 ? "|" : "");
        }

        out << "  leaf_" << leaf_index << " [label=\"{" << fields
            << "}\", shape=record];\n";
    }

    /**
     * @brief Checks if the high order bit of a value is set.
     * @param value The value to be checked.
     * @return True if the high order bit is set, false otherwise.
     */
    static constexpr bool isHighOrderBitSet(uint32_t value)
    {
        return value & (1u << 31);
    }

    /**
     * @brief Clears the high order bit of a value.
     * @param value The value to be modified.
     * @return The value with the high order bit cleared.
     */
    static constexpr uint32_t clearHighOrderBit(uint32_t value)
    {
        return value & ~(1u << 31);
    }

    /**
     * @brief Replaces the ".dot" extension in a filename with ".png".
     * @param filename The original filename.
     * @return The modified filename with the ".png" extension.
     */
    std::string replaceDotWithPng(const std::string &filename) const
    {
        std::string newFilename = filename;
        size_t pos = newFilename.rfind(".dot"); // Search for .dot in the filename
        if (pos != std::string::npos)
        {
            newFilename.replace(pos, 4, ".png"); // Replace .dot with .png
        }
        return newFilename;
    }

    constexpr static unsigned Terminal = 0xFFFFFFFF; ///< Marks that the bin is terminal.
    static constexpr unsigned Filler = 0xFFFFFFFE;   ///< Marks freed space in the vectors.

    const std::vector<uint32_t> &inner_nodes_; ///< Vector of inner node values.
    const std::vector<uint32_t> &leaf_nodes_;  ///< Vector of leaf node values.
    size_t num_bins_;                          ///< Number of bins for each node.
};