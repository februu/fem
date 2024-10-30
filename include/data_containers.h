#pragma once

const int NODES_PER_ELEMENT = 4;
const int NUMBER_OF_INTEGRATION_POINTS = 4;

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
    Jacobian jacobians[NUMBER_OF_INTEGRATION_POINTS];
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

    void print()
    {
        std::cout << "\n=== Nodes: ===\n";
        for (int i = 0; i < amountOfNodes; ++i)
            std::cout << "Node " << (i + 1) << ": ("
                      << nodes[i].x << ", "
                      << nodes[i].y << ")\n";

        std::cout << "\n=== Elements: ===\n";
        for (int i = 0; i < amountOfElements; ++i)
            std::cout << "Element " << (i + 1) << ": ("
                      << elements[i].nodeIds[0] + 1 << ", "
                      << elements[i].nodeIds[1] + 1 << ", "
                      << elements[i].nodeIds[2] + 1 << ", "
                      << elements[i].nodeIds[3] + 1 << ")\n";
    }
};

struct UniversalElement
{
    double dN_dKsi[NUMBER_OF_INTEGRATION_POINTS][4];
    double dN_dEta[NUMBER_OF_INTEGRATION_POINTS][4];

    void initialize();
};
