//
// Created by lucasalf on 15-05-2024.
//

#include "GraphLoader.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
/**
  * @brief Constructor for GraphLoader.
  *
  * @param g Pointer to the graph where the data will be loaded.
  */
GraphLoader::GraphLoader(Graph<int>* g) {
    this->g = g;
}

//------- Toy Graphs -------

/**
  * @brief Loads a toy graph from a CSV file.
  *
  * This function reads data from a CSV file containing information about toy graphs.
  * The data includes source, destination, and weight for each edge.
  *
  * @param path Path to the CSV file containing the toy graph data.
   */
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
        int source, dest;
        double d;

        getline(ss, line, ',');
        source = stoi(line);
        getline(ss, line, ',');
        dest = stoi(line);

        if (path != "../data/Toy-Graphs/tourism.csv"){

            getline(ss, line, '\r');
            d = stod(line);

            if (g->findVertex(source) == nullptr) {
                g->addVertex(source);
            }
            if (g->findVertex(dest) == nullptr) {
                g->addVertex(dest);
            }
            g->addBidirectionalEdge(source, dest, d);
        }
        else{
            string label1;
            string label2;

            getline(ss, line, ',');
            d = stod(line);
            getline(ss, label1, ',');
            getline(ss, label2, '\r');

            if (g->findVertex(source) == nullptr) {
                g->addVertex(source);
                g->findVertex(source)->setLabel(label1);
            }
            if (g->findVertex(dest) == nullptr) {
                g->addVertex(dest);
                g->findVertex(dest)->setLabel(label2);
            }
            g->addBidirectionalEdge(source, dest, d);
        }
    }
    file.close();
}

//------- Fully-Connected Graphs -------
/**
  * @brief Loads a fully connected graph from a CSV file.
  *
  * This function reads data from a CSV file containing information about fully connected graphs.
  * The data includes source, target, and weight for each edge.
  *
  * @param path Path to the directory containing the CSV file with the graph data.
  * @param size Size of the graph (number of vertices).
  */
void GraphLoader::loadFullyConnectedGraph(const std::string path, std::string size){
    ifstream file(path + "edges_" + size + ".csv");

    if (!file.is_open()){
        cout << "Error opening file" << endl;
        return;
    }

    string line;
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
            g->addVertex(source);
        }
        if (g->findVertex(target) == nullptr) {
            g->addVertex(target);
        }
        g->addBidirectionalEdge(source, target, d);
    }
}

//------- Real Graphs -------

/**
 * @brief Loads a real-world graph from CSV files.
 *
 * This function loads vertex and edge data from separate CSV files for real-world graphs.
 *
 * @param path Path to the directory containing the CSV files with vertex and edge data.
 */
void GraphLoader::loadRealGraph(const std::string& path){
    loadVertexReal(path + "/nodes.csv");
    loadEdgesReal(path + "/edges.csv");
}

/**
  * @brief Loads vertex data for a real-world graph from a CSV file.
  *
  * This function reads vertex data from a CSV file containing information about vertices.
  * The data includes vertex ID, longitude, and latitude.
  *
  * @param path Path to the CSV file containing the vertex data.
  */
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
    file.close();
}

/**
 * * @brief Loads edge data for a real-world graph from a CSV file.
 *
 * This function reads edge data from a CSV file containing information about edges.
 * The data includes source, target, and weight for each edge.
 *
 * @param path Path to the CSV file containing the edge data.
 */
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
    file.close();
}