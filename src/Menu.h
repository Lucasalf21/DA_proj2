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

class Menu {
    Graph<int>* g;

public:
    Menu(Graph<int>* g);
    void mainMenu();
    void choice1();
    void choice4();

    bool isGraphConnected();
    void dfs(int v, std::unordered_set<int>& visited);

    void TSPBacktracking(int curr, std::vector<int>& path, std::vector<int>& bestPath, double& bestDist, double currDist);
    double TSPNearestNeighbor(int startNode, std::vector<int>& path);
};


#endif //PROJ2_2324_MENU_H
