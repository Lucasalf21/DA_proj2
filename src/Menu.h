//
// Created by lucasalf on 15-05-2024.
//

#ifndef PROJ2_2324_MENU_H
#define PROJ2_2324_MENU_H

#include <iostream>
#include <chrono>
#include "data_structures/Graph.h"
#include "GraphLoader.h"

class Menu {
    Graph<int>* g;

public:
    Menu(Graph<int>* g);
    void mainMenu();
    void choice1();
    void choice3();

    void TSPBacktracking(int curr, std::vector<int>& path, std::vector<int>& bestPath, double& bestDist, double currDist);
    static double nearestNeighborTSP(Graph<int>* graph);
    static int nearestNeighbor(const Graph<int>* graph, const std::vector<bool>& visited, int current, double& minDistance, bool isreal);
    static double haversineDistance(const Vertex<int>& node1, const Vertex<int>& node2);
};


#endif //PROJ2_2324_MENU_H
