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
    {
        int gaussIndexKsi = i == 0 || i == 3 ? 0 : 1;
        int gaussIndexEta = i / 2;

        dN_dKsi[i][0] = -0.25 * (1 - gaussNodes[1][gaussIndexEta]);
        dN_dKsi[i][1] = 0.25 * (1 - gaussNodes[1][gaussIndexEta]);
        dN_dKsi[i][2] = 0.25 * (1 + gaussNodes[1][gaussIndexEta]);
        dN_dKsi[i][3] = -0.25 * (1 + gaussNodes[1][gaussIndexEta]);

        dN_dEta[i][0] = -0.25 * (1 - gaussNodes[1][gaussIndexKsi]);
        dN_dEta[i][1] = -0.25 * (1 + gaussNodes[1][gaussIndexKsi]);
        dN_dEta[i][2] = 0.25 * (1 + gaussNodes[1][gaussIndexKsi]);
        dN_dEta[i][3] = 0.25 * (1 - gaussNodes[1][gaussIndexKsi]);
    }
}