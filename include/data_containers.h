#pragma once

#include "quad.h"

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
    bool isBoundary = false;

    Node() : x(0.0), y(0.0) {}
    Node(double x, double y) : x(x), y(y) {}
};

struct Element
{
    int nodeIds[4];
    Jacobian jacobians[NUMBER_OF_INTEGRATION_POINTS_2D];
    double dN_dx[NUMBER_OF_INTEGRATION_POINTS_2D][4];
    double dN_dy[NUMBER_OF_INTEGRATION_POINTS_2D][4];
    double H[4][4] = {0};
    double C[4][4] = {0};
    double Hbc[4][4] = {0};
    double P[4] = {0};

    void printH();
    void printC();
    void printHbc();
    void printP();
};

struct Surface
{
    double N[NUMBER_OF_INTEGRATION_POINTS][4] = {0};
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

    // Copy constructor
    Grid(const Grid &other) : amountOfNodes(other.amountOfNodes), amountOfElements(other.amountOfElements)
    {
        nodes = new Node[amountOfNodes];
        elements = new Element[amountOfElements];

        for (int i = 0; i < amountOfNodes; i++)
            nodes[i] = other.nodes[i];

        for (int i = 0; i < amountOfElements; i++)
            elements[i] = other.elements[i];
    }

    Grid &operator=(const Grid &other)
    {
        if (this != &other)
        {
            delete[] nodes;
            delete[] elements;

            amountOfNodes = other.amountOfNodes;
            amountOfElements = other.amountOfElements;

            nodes = new Node[amountOfNodes];
            elements = new Element[amountOfElements];

            for (int i = 0; i < amountOfNodes; i++)
                nodes[i] = other.nodes[i];

            for (int i = 0; i < amountOfElements; i++)
                elements[i] = other.elements[i];
        }
        return *this;
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

    double N[NUMBER_OF_INTEGRATION_POINTS_2D][4];

    Surface surfaces[4];

    void initialize();
    void print();
};

struct Solution
{
    double **H, **C, *P, *T;
    int amountOfNodes;
    int simulationStep;

    Solution(int amountOfNodes, int initialTemp, int simulationStep) : amountOfNodes(amountOfNodes), simulationStep(simulationStep)
    {
        H = new double *[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            H[i] = new double[amountOfNodes]();

        C = new double *[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            C[i] = new double[amountOfNodes]();

        P = new double[amountOfNodes]();

        T = new double[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            T[i] = initialTemp;
    }

    Solution(const Solution &other) : amountOfNodes(other.amountOfNodes)
    {
        H = new double *[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            H[i] = new double[amountOfNodes];

        C = new double *[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            C[i] = new double[amountOfNodes];

        for (int i = 0; i < amountOfNodes; i++)
            for (int j = 0; j < amountOfNodes; j++)
                H[i][j] = other.H[i][j];

        for (int i = 0; i < amountOfNodes; i++)
            for (int j = 0; j < amountOfNodes; j++)
                C[i][j] = other.C[i][j];

        P = new double[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            P[i] = other.P[i];

        T = new double[amountOfNodes];
        for (int i = 0; i < amountOfNodes; i++)
            T[i] = other.T[i];
    }

    ~Solution()
    {
        for (int i = 0; i < amountOfNodes; i++)
            delete[] H[i];
        for (int i = 0; i < amountOfNodes; i++)
            delete[] C[i];
        delete[] H;
        delete[] C;
        delete[] P;
        delete[] T;
    }

    void printH();
    void printC();
    void printP();
    void printT();
    void printTMinMax();
    void solve();
};
