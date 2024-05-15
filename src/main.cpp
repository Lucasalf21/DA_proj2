#include <iostream>
#include "data_structures/Graph.h"
#include "GraphLoader.h"

using namespace std;

int main() {
    Graph<int> g;
    GraphLoader loader(&g);
    loader.loadToyGraph("../data/Toy-Graphs/tourism.csv");
    for (auto v : g.getVertexSet()) {
        cout << v->getInfo() << "-" << v->getLabel() << endl;
    }
    return 0;
}
