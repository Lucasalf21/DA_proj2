//
// Created by lucasalf on 15-05-2024.
//

#include "GraphLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

GraphLoader::GraphLoader(Graph<int>* g) {
    this->g = g;
}

void GraphLoader::loadToyGraph(const std::string& path){
    ifstream file(path);

    if (!file.is_open()){
        cout << "Error opening file" << endl;
        return;
    }

    string line;
    getline(file, line);

    while (getline(file, line)){
        stringstream ss(line);
        int v1, v2;
        double d;

        getline(ss, line, ',');
        v1 = stoi(line);
        getline(ss, line, ',');
        v2 = stoi(line);

        if (path != "../data/Toy-Graphs/tourism.csv"){

            getline(ss, line, '\r');
            d = stod(line);

            if (g->findVertex(v1) == nullptr) {
                g->addVertex(v1);
            }
            if (g->findVertex(v2) == nullptr) {
                g->addVertex(v2);
            }
            g->addEdge(v1, v2, d);
        }
        else{
            string label1;
            string label2;

            getline(ss, line, ',');
            d = stod(line);
            getline(ss, label1, ',');
            getline(ss, label2, '\r');

            if (g->findVertex(v1) == nullptr) {
                g->addVertex(v1);
                g->findVertex(v1)->setLabel(label1);
            }
            if (g->findVertex(v2) == nullptr) {
                g->addVertex(v2);
                g->findVertex(v2)->setLabel(label2);
            }
            g->addEdge(v1, v2, d);
        }
    }

    file.close();
}

void GraphLoader::loadVertexReal(std::string path) {
    ifstream file(path);

    if (!file.is_open()){
        cout << "Error opening file" << endl;
        return;
    }

    string line;
    getline(file, line);

    while (getline(file, line)){
        stringstream ss(line);
        int id;
        double lon, lat;

        getline(ss, line, ',');
        id = stoi(line);
        getline(ss, line, ',');
        lon = stod(line);
        getline(ss, line, '\r');
        lat = stod(line);

        if (g->findVertex(id) == nullptr) {
            g->addVertex(id);
            g->findVertex(id)->setLongitude(lon);
            g->findVertex(id)->setLatitude(lat);
        }
    }
}

void GraphLoader::loadEdgesReal(std::string path) {
    ifstream file(path);

    if (!file.is_open()){
        cout << "Error opening file" << endl;
        return;
    }

    string line;
    getline(file, line);

    while (getline(file, line)){
        stringstream ss(line);
        int source, target;
        double d;

        getline(ss, line, ',');
        source = stoi(line);
        getline(ss, line, ',');
        target = stoi(line);
        getline(ss, line, '\r');
        d = stod(line);

        if (g->findVertex(source) == nullptr) {
            cout << "Error loading edge: source vertex not found" << endl;
            return;
        }
        if (g->findVertex(target) == nullptr) {
            cout << "Error loading edge: target vertex not found" << endl;
            return;
        }
        g->addEdge(source, target, d);
    }
}