//
// Created by lucasalf on 15-05-2024.
//

#include "Menu.h"
#include <algorithm>
#include <stack>
#include <cmath>

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
                choice4();
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
    if (option == 'n') {
        return;
    }
    if (g->getVertexSet().size() > 11) {
        cout << "This algorithm will take a REALLY long time for this graph. Do you want to continue? (y/n)" << endl;
        cin >> option;
        cout << endl;
        if (option == 'n') {
            return;
        }
    }
    vector<int> path = {0};
    vector<int> bestPath;
    double minDist = numeric_limits<double>::infinity();
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

void Menu::choice4() {
    int startNode;
    cout << "Enter the starting node for the TSP tour: ";
    cin >> startNode;
    cout << endl;

    unordered_set<int> visited;
    if (!isGraphConnected()) {
        cout << "No valid TSP tour exists as the graph is not fully connected from the start node." << endl;
        return;
    }

    vector<int> path;
    auto start = chrono::high_resolution_clock::now();
    double totalDist = TSPNearestNeighbor(startNode, path);
    auto end = chrono::high_resolution_clock::now();

    if (totalDist == numeric_limits<double>::infinity()) {
        cout << "No valid TSP tour exists." << endl;
    } else {
        cout << "Best Path: " << endl;
        for (int i = 0; i < path.size() - 1; i++) {
            cout << path[i] << " -> ";
        }
        cout << path.back() << endl;
        cout << "Distance: " << totalDist << endl;
        cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl << endl;
    }
}

double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371.0;

    lat1 = lat1 * M_PI / 180.0;
    lon1 = lon1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;
    lon2 = lon2 * M_PI / 180.0;

    // Differences
    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    // Haversine formula
    double a = pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return R * c;
}

bool Menu::isGraphConnected() {
    unordered_set<int> visited;
    vector<int> components;
    for (auto v : g->getVertexSet()) {
        if (visited.find(v.second->getInfo()) == visited.end()) {
            dfs(v.second->getInfo(), visited);
            components.push_back(v.second->getInfo());
        }
    }
    for (int i = 0; i < components.size() - 1; i++) {
        double lat1 = g->findVertex(components[i])->getLatitude();
        double lon1 = g->findVertex(components[i])->getLongitude();
        double lat2 = g->findVertex(components[i+1])->getLatitude();
        double lon2 = g->findVertex(components[i+1])->getLongitude();
        double weight = haversine(lat1, lon1, lat2, lon2);
        g->addEdge(components[i], components[i+1], weight);
    }
    return visited.size() == g->getVertexSet().size();
}

void Menu::dfs(int v, unordered_set<int>& visited) {
    visited.insert(v);
    for (auto edge : g->findVertex(v)->getAdj()) {
        int neighbor = edge.second->getDest()->getInfo();
        if (visited.find(neighbor) == visited.end()) {
            dfs(neighbor, visited);
        }
    }
}

double Menu::TSPNearestNeighbor(int startNode, vector<int>& path) {
    unordered_set<int> visited;
    std::stack<int> pathStack;
    int currentNode = startNode;
    double totalDist = 0;
    path.push_back(currentNode);
    pathStack.push(currentNode);
    visited.insert(currentNode);

    while (!pathStack.empty()) {
        currentNode = pathStack.top();
        auto currVertex = g->findVertex(currentNode);
        double minDist = numeric_limits<double>::infinity();
        int nextNode = -1;

        for (auto e : currVertex->getAdj()) {
            int neighbor = e.second->getDest()->getInfo();
            if (visited.find(neighbor) == visited.end() && e.second->getWeight() < minDist) {
                minDist = e.second->getWeight();
                nextNode = neighbor;
            }
        }

        if (nextNode == -1) {
            pathStack.pop();
        } else {
            currentNode = nextNode;
            pathStack.push(currentNode);
            path.push_back(currentNode);
            visited.insert(currentNode);
            totalDist += minDist;
        }
    }

    auto lastEdge = g->findEdge(currentNode, startNode);
    if (lastEdge != nullptr) {
        totalDist += lastEdge->getWeight();
    } else {
        path.pop_back();
    }

    return totalDist;
}

void Menu::TSPBacktracking(int curr, vector<int>& path, vector<int>& bestPath, double& bestDist, double currDist) {
    auto currVertex = g->findVertex(curr);
    currVertex->setVisited(true);
    if (path.size() == g->getVertexSet().size()) {
        double totalDist = currDist + g->findEdge(curr, path[0])->getWeight();
        if (totalDist < bestDist) {
            bestDist = totalDist;
            bestPath = path;
        }
        currVertex->setVisited(false);
        return;
    }
    for (auto e : currVertex->getAdj()) {
        auto w = e.second->getDest();
        if (!w->isVisited()) {
            double newDist = currDist + e.second->getWeight();
            if (newDist < bestDist) {
                path.push_back(w->getInfo());
                TSPBacktracking(w->getInfo(), path, bestPath, bestDist, newDist);
                path.pop_back();
            }
        }
    }
    currVertex->setVisited(false);
}
