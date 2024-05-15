//
// Created by lucasalf on 15-05-2024.
//

#include "Menu.h"

using namespace std;

Menu::Menu(Graph<int>* g) {
    this->g = g;
    mainMenu();
}

void Menu::mainMenu() {
    int option;
    cout << "|---------------------------------------------|" << endl
         << "|     Ocean Shipping and Urban Deliveries     |" << endl
         << "|---------------------------------------------|" << endl << endl

         << "1 - Calculate TSP using a backtracking approach (4.1)." << endl
         << "2 - Triangular Approximation Heuristic (4.2)." << endl
         << "3 - Other Heuristics (4.3)." << endl
         << "4 - TSP in the Real World (4.4)." << endl << endl

         << "0 - Exit" << endl << endl;

    cout << "Option: ";
    cin >> option;

    switch(option) {
        case 1:
            TSPBacktracking();
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

void Menu::TSPBacktracking(){

}