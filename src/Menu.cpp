//
// Created by lucasalf on 15-05-2024.
//

#include "Menu.h"
#include <algorithm>
#include <stack>
#include <cmath>
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
    std::chrono::duration<double> duration = end - start;

    cout << "Best Path: " << endl;
    for (auto v : bestPath) {
        cout << v << " -> ";
    }
    cout << bestPath[0] << endl;
    cout << "Distance: " << minDist << endl;
    cout << "Time: " << duration.count() << " seconds\n";
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
    std::chrono::duration<double> duration = end - start;

    if (totalDist == numeric_limits<double>::infinity()) {
        cout << "No valid TSP tour exists." << endl;
    } else {
        cout << "Best Path: " << endl;
        for (int i = 0; i < path.size() - 1; i++) {
            cout << path[i] << " -> ";
        }
        cout << path.back() << endl;
        cout << "Distance: " << totalDist << endl;
        cout << "Time: " << duration.count() << " seconds" << endl << endl;
    }
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
 * @brief Checks whether the graph is fully connected.
 *
 * @return True if the graph is fully connected, false otherwise.
 */
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
        double weight = haversineDistance(*g->findVertex(components[i]), *g->findVertex(components[i+1]));
        g->addEdge(components[i], components[i+1], weight);
    }
    return visited.size() == g->getVertexSet().size();
}

/**
 * @brief Performs depth-first search (DFS) traversal on the graph to check for connectivity.
 *
 * @param v The starting vertex for DFS traversal.
 * @param visited A set to keep track of visited vertices.
 */
void Menu::dfs(int v, unordered_set<int>& visited) {
    visited.insert(v);
    for (auto edge : g->findVertex(v)->getAdj()) {
        int neighbor = edge.second->getDest()->getInfo();
        if (visited.find(neighbor) == visited.end()) {
            dfs(neighbor, visited);
        }
    }
}

/**
 * @brief Finds the best TSP tour starting from a specified node using the nearest neighbor heuristic.
 *
 * @param startNode The starting node for the TSP tour.
 * @param path[out] The vector to store the TSP tour path.
 * @return The total distance of the TSP tour.
 */
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
    bool isRealGraph = graph->getVertexSet().begin()->second->getLongitude() != 0;
    bool found = false;

    Graph<int> mstGraph;
    for (auto vertex : graph->getVertexSet()) {
        mstGraph.addVertex(vertex.second->getInfo());
    }
    for (auto edge : mstEdges) {
        mstGraph.addEdge(edge->getOrig()->getInfo(), edge->getDest()->getInfo(), edge->getWeight());
        mstGraph.addEdge(edge->getDest()->getInfo(), edge->getOrig()->getInfo(), edge->getWeight());
    }

    dfsMST(mstGraph.findVertex(0), visited, tour);

    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        Vertex<int>* u = graph->findVertex(tour[i]);
        Vertex<int>* v = graph->findVertex(tour[i + 1]);
        for (auto e : u->getAdj()) {
            if (e.second->getDest() == v) {
                totalDistance += e.second->getWeight();
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
        if (e.second->getDest() == graph->findVertex(tour.front())) {
            totalDistance += e.second->getWeight();
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

    for (auto edge : vertex->getAdj()) {
        Vertex<int>* dest = edge.second->getDest();
        if (visited.find(dest->getInfo()) == visited.end()) {
            dfsMST(dest, visited, tour);
            tour.push_back(vertex->getInfo());
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

    for (auto vertex : vertices) {
        vertex.second->setDist(std::numeric_limits<double>::infinity());
        vertex.second->setPath(nullptr);
    }

    Vertex<int>* startVertex = vertices[0];
    startVertex->setDist(0);
    pq.insert(startVertex);

    while (!pq.empty()) {
        Vertex<int>* u = pq.extractMin();
        if (inMST.find(u->getInfo()) != inMST.end()) continue;
        inMST.insert(u->getInfo());

        for (auto edge : u->getAdj()) {
            Vertex<int>* v = edge.second->getDest();
            double weight = edge.second->getWeight();

            if (inMST.find(v->getInfo()) == inMST.end() && weight < v->getDist()) {
                v->setDist(weight);
                v->setPath(edge.second);
                pq.insert(v);
            }
        }
    }

    // Collect MST edges
    for (auto vertex : vertices) {
        if (vertex.second->getPath() != nullptr) {
            mstEdges.push_back(vertex.second->getPath());
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
    bool isRealGraph = graph->getVertexSet().begin()->second->getLongitude() != 0;
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
        if (e.second->getDest() == graph->findVertex(0)) {
            totalDistance += e.second->getWeight();
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

    for (auto vertex : graph->getVertexSet()) {
        if (!visited[vertex.second->getInfo()]) {
            double distance;
            for (auto e : currentVertex->getAdj()) {
                if (e.second->getDest() == vertex.second) {
                    distance = e.second->getWeight();
                    found = true;
                    break;
                }
            }
            if (isRealGraph && !found) {
                distance = haversineDistance(*currentVertex, *vertex.second);
                found = false;
            }

            if (distance < minDistance) {
                minDistance = distance;
                nearest = vertex.second->getInfo();
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