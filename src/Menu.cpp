//
// Created by lucasalf on 15-05-2024.
//

#include "Menu.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>

using namespace std;

/**
 * @brief Constructs a new Menu object.
 *
 * @param g Pointer to the graph.
 */
Menu::Menu(Graph<int>* g) {
    this->g = g;
    mainMenu();
}

/**
 * @brief Displays the main menu and handles user input.
 */
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
                choice2();
                break;
            case 3:
                // Other Heuristics
                choice3();
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

/**
 * @brief Executes the backtracking approach for TSP.
 *
 * Time Complexity: O(n!), where n is the number of vertices.
 */
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
    std::chrono::duration<double> duration = end - start;

    cout << "Best Path: " << endl;
    for (auto v : bestPath) {
        cout << v << " -> ";
    }
    cout << bestPath[0] << endl;
    cout << "Distance: " << minDist << endl;
    cout << "Time: " << duration.count() << " seconds\n";
}

/**
 * @brief Executes the triangular approximation heuristic for TSP.
 *
 * Time Complexity: O(n^2), where n is the number of vertices.
 */
void Menu::choice2(){
    auto start = std::chrono::high_resolution_clock::now();
    double totalDistance = triangularApproximationTSP(this->g);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Distance: " << totalDistance << std::endl;
    std::cout << "Time: " << duration.count() << " seconds\n";
}

/**
 * @brief Executes the nearest neighbor heuristic for TSP.
 *
 * Time Complexity: O(n^2), where n is the number of vertices.
 */
void Menu::choice3() {
    auto start = std::chrono::high_resolution_clock::now();
    double totalDistance = nearestNeighborTSP(this->g);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Distance: " << totalDistance << std::endl;
    std::cout << "Time: " << duration.count() << " seconds\n";
}

/**
 * @brief Recursive function for backtracking TSP.
 *
 * @param curr The current vertex.
 * @param path The current path.
 * @param bestPath The best path found.
 * @param bestDist The best distance found.
 * @param currDist The current distance.
 *
 * Time Complexity: O(n!), where n is the number of vertices.
 */
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

/**
 * @brief Executes the triangular approximation heuristic for TSP.
 *
 * @param graph Pointer to the graph.
 * @return The total distance of the tour.
 *
 * Time Complexity: O(n^2), where n is the number of vertices.
 */
double Menu::triangularApproximationTSP(Graph<int>* graph) {
    std::vector<Edge<int>*> mstEdges = primMST(graph);
    std::unordered_set<int> visited;
    std::vector<int> tour;
    bool isRealGraph = graph->getVertexSet().front()->getLongitude() != 0;
    bool found = false;

    // Reconstruct the graph with MST edges only
    Graph<int> mstGraph;
    for (auto vertex : graph->getVertexSet()) {
        mstGraph.addVertex(vertex->getInfo());
    }
    for (auto edge : mstEdges) {
        mstGraph.addEdge(edge->getOrig()->getInfo(), edge->getDest()->getInfo(), edge->getWeight());
        mstGraph.addEdge(edge->getDest()->getInfo(), edge->getOrig()->getInfo(), edge->getWeight());
    }

    // Perform DFS on the MST to get the TSP tour
    dfsMST(mstGraph.findVertex(0), visited, tour);

    // Calculate the total distance of the tour
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        Vertex<int>* u = graph->findVertex(tour[i]);
        Vertex<int>* v = graph->findVertex(tour[i + 1]);
        for (auto e : u->getAdj()) {
            if (e->getDest() == v) {
                totalDistance += e->getWeight();
                found = true;
                break;
            }
        }
        if (isRealGraph && !found) {
            totalDistance += haversineDistance(*u, *v);
        }
        found = false;
    }

    for (auto e : graph->findVertex(tour.back())->getAdj()) {
        if (e->getDest() == graph->findVertex(tour.front())) {
            totalDistance += e->getWeight();
            found = true;
        }
    }

    if (isRealGraph && !found) {
        totalDistance += haversineDistance(*graph->findVertex(tour.back()), *graph->findVertex(tour.front()));
    }

    return totalDistance;
}

/**
 * @brief Depth-first search on the MST to get the TSP tour.
 *
 * @param vertex The current vertex.
 * @param visited Set of visited vertices.
 * @param tour The TSP tour.
 *
 * Time Complexity: O(n), where n is the number of vertices.
 */
void Menu::dfsMST(Vertex<int>* vertex, std::unordered_set<int>& visited, std::vector<int>& tour) {
    visited.insert(vertex->getInfo());
    tour.push_back(vertex->getInfo());

    for (Edge<int>* edge : vertex->getAdj()) {
        Vertex<int>* dest = edge->getDest();
        if (visited.find(dest->getInfo()) == visited.end()) {
            dfsMST(dest, visited, tour);
            tour.push_back(vertex->getInfo()); // Return to the current vertex
        }
    }
}

/**
 * @brief Computes the minimum spanning tree using Prim's algorithm.
 *
 * @param graph Pointer to the graph.
 * @return A vector of edges in the MST.
 *
 * Time Complexity: O(n^2), where n is the number of vertices.
 */
std::vector<Edge<int>*> Menu::primMST(Graph<int>* graph) {
    std::vector<Edge<int>*> mstEdges;
    std::unordered_set<int> inMST;
    MutablePriorityQueue<Vertex<int>> pq;

    auto vertices = graph->getVertexSet();
    vertices.front()->setDist(0);
    pq.insert(vertices.front());

    while (!pq.empty()) {
        Vertex<int>* u = pq.extractMin();
        if (inMST.find(u->getInfo()) != inMST.end()) continue;
        inMST.insert(u->getInfo());

        for (Edge<int>* edge : u->getAdj()) {
            Vertex<int>* v = edge->getDest();
            double weight = edge->getWeight();

            if (inMST.find(v->getInfo()) == inMST.end() && weight < v->getDist()) {
                v->setDist(weight);
                v->setPath(edge);
                pq.insert(v);
            }
        }
    }

    for (Vertex<int>* vertex : vertices) {
        if (vertex->getPath() != nullptr) {
            mstEdges.push_back(vertex->getPath());
        }
    }

    return mstEdges;
}

/**
 * @brief Executes the nearest neighbor heuristic for TSP.
 *
 * @param graph Pointer to the graph.
 * @return The total distance of the tour.
 *
 * Time Complexity: O(n^2), where n is the number of vertices.
 */
double Menu::nearestNeighborTSP(Graph<int>* graph) {
    double totalDistance = 0.0;
    std::vector<bool> visited(graph->getNumVertex(), false);
    int current = 0;
    visited[current] = true;
    bool isRealGraph = graph->getVertexSet().front()->getLongitude() != 0;
    bool found = false;

    for (size_t i = 1; i < graph->getNumVertex(); ++i) {
        double minDistance;
        int next = nearestNeighbor(graph, visited, current, minDistance, isRealGraph);
        if (next == -1) break;
        totalDistance += minDistance;
        visited[next] = true;
        current = next;
    }

    for (auto e : graph->findVertex(current)->getAdj()) {
        if (e->getDest() == graph->findVertex(0)) {
            totalDistance += e->getWeight();
            found = true;
            break;
        }
    }

    // Return to the starting node
    if (isRealGraph && !found) {
        totalDistance += haversineDistance(*graph->findVertex(current), *graph->findVertex(0));
    }

    return totalDistance;
}

/**
 * @brief Finds the nearest neighbor for the nearest neighbor TSP heuristic.
 *
 * @param graph Pointer to the graph.
 * @param visited Vector of visited vertices.
 * @param current The current vertex.
 * @param minDistance The minimum distance to the next vertex.
 * @param isRealGraph Indicates if the graph has real-world coordinates.
 * @return The index of the nearest neighbor.
 *
 * Time Complexity: O(n), where n is the number of vertices.
 */
int Menu::nearestNeighbor(const Graph<int>* graph, const std::vector<bool>& visited, int current, double& minDistance, bool isRealGraph) {
    minDistance = std::numeric_limits<double>::max();
    int nearest = -1;
    bool found = false;
    Vertex<int>* currentVertex = graph->findVertex(current);

    for (Vertex<int>* vertex : graph->getVertexSet()) {
        if (!visited[vertex->getInfo()]) {
            double distance;
            for (auto e : currentVertex->getAdj()) {
                if (e->getDest() == vertex) {
                    distance = e->getWeight();
                    found = true;
                    break;
                }
            }
            if (isRealGraph && !found) {
                distance = haversineDistance(*currentVertex, *vertex);
                found = false;
            }

            if (distance < minDistance) {
                minDistance = distance;
                nearest = vertex->getInfo();
            }
        }
    }
    return nearest;
}

/**
 * @brief Calculates the Haversine distance between two vertices.
 *
 * @param node1 The first vertex.
 * @param node2 The second vertex.
 * @return The Haversine distance.
 *
 * Time Complexity: O(1)
 */
double Menu::haversineDistance(const Vertex<int>& node1, const Vertex<int>& node2) {
    const double R = 6371; // Earth radius in kilometers
    double lat1 = node1.getLatitude() * M_PI / 180.0;
    double lon1 = node1.getLongitude() * M_PI / 180.0;
    double lat2 = node2.getLatitude() * M_PI / 180.0;
    double lon2 = node2.getLongitude() * M_PI / 180.0;
    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}