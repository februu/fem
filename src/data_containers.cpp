#include <iostream>
#include "../include/data_containers.h"
#include "../include/quad.h"

void Jacobian::inverse()
{
    calculateDetJ();
    invJ[0][0] = J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] = J[0][0] / detJ;
}

double Jacobian::calculateDetJ()
{
    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    return detJ;
}

void Jacobian::printJ()
{
    std::cout << "| " << J[0][0] << " " << J[0][1] << " |\n";
    std::cout << "| " << J[1][0] << " " << J[1][1] << " |\n";
}

void Jacobian::printInvJ()
{
    std::cout << "| " << invJ[0][0] << " " << invJ[0][1] << " |\n";
    std::cout << "| " << invJ[1][0] << " " << invJ[1][1] << " |\n";
}

void Grid::print()
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

void UniversalElement::initialize()
{
    for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
        for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
        {
            int currentIndex = i * NUMBER_OF_INTEGRATION_POINTS + j;

            dN_dKsi[currentIndex][0] = -0.25 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            dN_dKsi[currentIndex][1] = 0.25 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            dN_dKsi[currentIndex][2] = 0.25 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            dN_dKsi[currentIndex][3] = -0.25 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);

            dN_dEta[currentIndex][0] = -0.25 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][i]);
            dN_dEta[currentIndex][1] = -0.25 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][i]);
            dN_dEta[currentIndex][2] = 0.25 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][i]);
            dN_dEta[currentIndex][3] = 0.25 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][i]);
        }
}

void Solution::printH()
{
    std::cout << "\n=== H (Global): ===\n";
    for (int i = 0; i < amountOfNodes; i++)
    {
        for (int j = 0; j < amountOfNodes; j++)
            std::cout << H[i][j] << " ";
        std::cout << "\n";
    }
}