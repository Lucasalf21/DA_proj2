//
// Created by lucasalf on 15-05-2024.
//

#ifndef PROJ2_2324_MENU_H
#define PROJ2_2324_MENU_H

#include <iostream>
#include <chrono>
#include <unordered_set>
#include "data_structures/Graph.h"
#include "GraphLoader.h"
/**
 * @class Menu
 * @brief This class provides various methods to solve the Traveling Salesman Problem (TSP) using different approaches.
 */
class Menu {
    Graph<int>* g;///< Pointer to the graph.

public:
    Menu(Graph<int>* g);
    void mainMenu();
    void choice1();
    void choice2();
    void choice3();

    void TSPBacktracking(int curr, std::vector<int>& path, std::vector<int>& bestPath, double& bestDist, double currDist);
    double triangularApproximationTSP(Graph<int>* graph);
    std::vector<Edge<int>*> primMST(Graph<int>* graph);
    void dfsMST(Vertex<int>* vertex, std::unordered_set<int>& visited, std::vector<int>& tour);
    static double nearestNeighborTSP(Graph<int>* graph);
    static int nearestNeighbor(const Graph<int>* graph, const std::vector<bool>& visited, int current, double& minDistance, bool isreal);
    static double haversineDistance(const Vertex<int>& node1, const Vertex<int>& node2);
};


#endif //PROJ2_2324_MENU_H
