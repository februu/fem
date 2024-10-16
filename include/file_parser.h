#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

const int NODES_PER_ELEMENT = 4;

struct Node
{
    double x, y;
};

struct Element
{
    int nodeIds[NODES_PER_ELEMENT];
};

struct Grid
{
    int amountOfNodes;
    int amountOfElements;
    Node *nodes;
    Element *elements;
};

struct GlobalData
{
    int simulationTime;
    int simulationStep;
    int conductivity;
    int alfa;
    int tot;
    int initialTemperature;
    int density;
    int specificHeat;
    Grid grid;

    void parseFile(const std::string &fileName);
    void print();
    ~GlobalData();
};
