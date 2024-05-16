//
// Created by lucasalf on 15-05-2024.
//

#ifndef PROJ2_2324_GRAPHLOADER_H
#define PROJ2_2324_GRAPHLOADER_H

#include "data_structures/Graph.h"

class GraphLoader {
private:
    Graph<int>* g;

public:
    GraphLoader(Graph<int>* g);
    void loadToyGraph(const std::string& path);
    void loadRealGraph(const std::string& path);
    void loadFullyConnectedGraph(const std::string path, std::string size);

    void loadVertexReal(std::string path);
    void loadEdgesReal(std::string path);
};


#endif //PROJ2_2324_GRAPHLOADER_H
