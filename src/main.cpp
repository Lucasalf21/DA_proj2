#include <iostream>
#include <set>
#include "data_structures/Graph.h"
#include "GraphLoader.h"
#include "Menu.h"

using namespace std;

void loadToyGraph(Graph<int>* g, GraphLoader loader) {
    set<int> options = {0, 1, 2, 3};
    cout << "Please select which toy graph you want to load:" << endl
         << "1 - Shipping" << endl
         << "2 - Stadiums" << endl
         << "3 - Tourism" << endl << endl

         << "0 - Exit" << endl << endl

         << "Option: ";

    int option;
    cin >> option;

    cout << endl;

    if (options.find(option) == options.end()) {
        cout << "Invalid option!" << endl;
        exit(0);
    }

    switch (option) {
        case 1:
            loader.loadToyGraph("../data/Toy-Graphs/shipping.csv");
            break;
        case 2:
            loader.loadToyGraph("../data/Toy-Graphs/stadiums.csv");
            break;
        case 3:
            loader.loadToyGraph("../data/Toy-Graphs/tourism.csv");
            break;
        case 0:
            exit(0);
    }
}

void loadFullyConnectedGraph(Graph<int>* g, GraphLoader loader) {
    std::string size;
    cout << "Please select how many vertices you want in the graph (25, 50, 75, 100, 200, 300, 400, 500, 600, 700, 800, 900) or 0 to exit: ";
    cin >> size;

    if (size == "0") {
        exit(0);
    }

    loader.loadFullyConnectedGraph("../data/Extra_Fully_Connected_Graphs/", size);
}

void loadRealGraph(Graph<int>* g, GraphLoader loader) {
    set<int> options = {0, 1, 2, 3};
    cout << "Please select which real graph you want to load:" << endl
         << "1 - Graph 1" << endl
         << "2 - Graph 2" << endl
         << "3 - Graph 3" << endl << endl

         << "0 - Exit" << endl << endl

         << "Option: ";

    int option;
    cin >> option;

    cout << endl;

    if (options.find(option) == options.end()) {
        cout << "Invalid option!" << endl;
        exit(0);
    }

    switch (option) {
        case 1:
            loader.loadRealGraph("../data/Real-world Graphs/graph1");
            break;
        case 2:
            loader.loadRealGraph("../data/Real-world Graphs/graph2");
            break;
        case 3:
            loader.loadRealGraph("../data/Real-world Graphs/graph3");
            break;
        case 0:
            exit(0);
    }
}

int main() {
    Graph<int>* g = new Graph<int>();
    GraphLoader loader(g);
    set<int> options = {0, 1, 2, 3};

    cout << "---------------------------------------------" << endl
         << "|               Graph Loader                |" << endl
         << "---------------------------------------------" << endl << endl

         << "Hello! Please select which graph you want to load:" << endl
         << "1 - Toy Graph" << endl
         << "2 - Fully-Connected graphs" << endl
         << "3 - Real Graph" << endl << endl

         << "0 - Exit" << endl << endl

         << "Option: ";

    int option;
    cin >> option;

    cout << endl;

    if (options.find(option) == options.end()) {
        cout << "Invalid option!" << endl;
        return 0;
    }

    switch (option) {
        case 1:
            loadToyGraph(g, loader);
            break;
        case 2:
            loadFullyConnectedGraph(g, loader);
            break;
        case 3:
            loadRealGraph(g, loader);
            break;
        case 0:
            exit(0);
    }

    cout << "Graph loaded!" << endl << endl;

    int edgeCount = 0;
    auto vertexSet = g->getVertexSet();
    for (auto &pair : vertexSet) {
        Vertex<int> *v = pair.second;
        for (auto e : v->getAdj()) {
            edgeCount++;
        }
    }

    Menu menu(g);

    return 0;
}
