//
// Created by lucasalf on 15-05-2024.
//

#include "Menu.h"
#include <algorithm>

using namespace std;

Menu::Menu(Graph<int>* g) {
    this->g = g;
    mainMenu();
}

void Menu::mainMenu() {
    int option;
    while (true) {
        cout << "|---------------------------------------------|" << endl
             << "|     Ocean Shipping and Urban Deliveries     |" << endl
             << "|---------------------------------------------|" << endl << endl

             << "1 - TSP using a backtracking approach (4.1)." << endl
             << "2 - Triangular Approximation Heuristic (4.2)." << endl
             << "3 - Other Heuristics (4.3)." << endl
             << "4 - TSP in the Real World (4.4)." << endl << endl

             << "0 - Exit" << endl << endl;

        cout << "Option: ";
        cin >> option;
        cout << endl;

        switch (option) {
            case 1:
                choice1();
                break;
            case 2:
                // Triangular Approximation Heuristic
                break;
            case 3:
                // Other Heuristics
                break;
            case 4:
                // TSP in the Real World
                break;
            case 0:
                exit(0);
            default:
                cout << "Invalid option" << endl;
                break;
        }
    }
}

void Menu::choice1() {
    cout << "WARNING - This algorithm will not work with the shipping toy graph. Do you wish to continue? (y/n)" << endl;
    char option;
    cin >> option;
    cout << endl;
    if (option == 'n'){
        return;
    }
    if (g->getVertexSet().size() > 11){
        cout << "This algorithm will take a REALLY long time for this graph. Do you want to continue? (y/n)" << endl;
        cin >> option;
        cout << endl;
        if (option == 'n'){
            return;
        }
    }

    vector<int> path = {0};
    vector<int> bestPath;
    double minDist = INF;

    auto start = chrono::high_resolution_clock::now();
    TSPBacktracking(0, path, bestPath, minDist, 0);
    auto end = chrono::high_resolution_clock::now();

    cout << "Best Path: " << endl;
    for (auto v : bestPath) {
        cout << v << " -> ";
    }
    cout << bestPath[0] << endl;
    cout << "Distance: " << minDist << endl;
    cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl << endl;
}

void Menu::TSPBacktracking(int curr, vector<int>& path, vector<int>& bestPath, double& bestDist, double currDist) {
    auto currVertex = g->findVertex(curr);
    currVertex->setVisited(true);

    if (path.size() == g->getVertexSet().size()) {
        double totalDist = currDist + g->findEdge(curr, 0)->getWeight();
        if (totalDist < bestDist) {
            bestDist = totalDist;
            bestPath = path;
        }
        return;
    }

    for (auto e : currVertex->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            w->setVisited(true);
            double newDist = currDist + e->getWeight();
            path.push_back(e->getDest()->getInfo());
            if (newDist < bestDist) {
                TSPBacktracking(w->getInfo(), path, bestPath, bestDist, newDist);
            }
            path.pop_back();
            w->setVisited(false);
        }
    }
}

