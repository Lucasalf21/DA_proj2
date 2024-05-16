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

    void TSPBacktracking(int curr, std::vector<int>& path, std::vector<int>& bestPath, double& bestDist, double currDist);
};


#endif //PROJ2_2324_MENU_H
