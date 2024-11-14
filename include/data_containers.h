#pragma once

#include <cmath>

const int NODES_PER_ELEMENT = 4;
const int NUMBER_OF_INTEGRATION_POINTS = 2;
const int NUMBER_OF_INTEGRATION_POINTS_2D = NUMBER_OF_INTEGRATION_POINTS * NUMBER_OF_INTEGRATION_POINTS;

struct Jacobian
{
    double J[2][2];
    double invJ[2][2];
    double detJ;

    void inverse();
    double calculateDetJ();
    void printJ();
    void printInvJ();
};

struct Node
{
    double x, y;

    Node() : x(0.0), y(0.0) {}
    Node(double x, double y) : x(x), y(y) {}
};

struct Element
{
    int nodeIds[NODES_PER_ELEMENT];
    Jacobian jacobians[NUMBER_OF_INTEGRATION_POINTS_2D];
    double dN_dx[NUMBER_OF_INTEGRATION_POINTS_2D][4];
    double dN_dy[NUMBER_OF_INTEGRATION_POINTS_2D][4];
    double H[4][4] = {0};
};

struct Grid
{
    int amountOfNodes;
    int amountOfElements;
    Node *nodes;
    Element *elements;

    Grid() : amountOfNodes(0), amountOfElements(0), nodes(nullptr), elements(nullptr) {}
    Grid(int amountOfNodes, int amountOfElements) : amountOfNodes(amountOfNodes), amountOfElements(amountOfElements)
    {
        nodes = new Node[amountOfNodes];
        elements = new Element[amountOfElements];
    }

    ~Grid()
    {
        delete[] nodes;
        delete[] elements;
    }

    void print();
};

struct UniversalElement
{
    double dN_dKsi[NUMBER_OF_INTEGRATION_POINTS_2D][4];
    double dN_dEta[NUMBER_OF_INTEGRATION_POINTS_2D][4];

    void initialize();
};
