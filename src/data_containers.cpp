#include <iostream>
#include "../include/data_containers.h"
#include "../include/quad.h"
#include "../include/gauss_elimination.h"

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
                  << nodes[i].y << "), isBoundary=" << nodes[i].isBoundary << "\n";

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
    // Calculate dN_dKsi and dN_dEta
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

    // Calculate N in integration points for each surface
    // 0 - bottom, 1 - right, 2 - top, 3 - left
    // Some N functions were skipped because they are multiplied by 0
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
        {
            surfaces[0].N[j][0] = 0.5 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            surfaces[0].N[j][1] = 0.5 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);

            surfaces[1].N[j][1] = 0.5 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            surfaces[1].N[j][2] = 0.5 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);

            surfaces[2].N[j][2] = 0.5 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            surfaces[2].N[j][3] = 0.5 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);

            surfaces[3].N[j][0] = 0.5 * (1 - gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
            surfaces[3].N[j][3] = 0.5 * (1 + gaussNodes[NUMBER_OF_INTEGRATION_POINTS - 1][j]);
        }
}

void UniversalElement::print()
{
    std::cout << "\n=== Universal Element: ===\n";
    std::cout << "dN_dKsi:\n";
    for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS_2D; i++)
    {
        for (int j = 0; j < 4; j++)
            std::cout << dN_dKsi[i][j] << " ";
        std::cout << "\n";
    }

    std::cout << "dN_dEta:\n";
    for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS_2D; i++)
    {
        for (int j = 0; j < 4; j++)
            std::cout << dN_dEta[i][j] << " ";
        std::cout << "\n";
    }

    std::cout << "Surfaces:\n";
    for (int i = 0; i < 4; i++)
    {
        std::cout << "Surface " << i << ":\n";
        for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
        {
            for (int k = 0; k < 4; k++)
                std::cout << surfaces[i].N[j][k] << " ";
            std::cout << "\n";
        }
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

void Solution::printP()
{
    std::cout << "\n=== P (Global): ===\n";
    for (int i = 0; i < amountOfNodes; i++)
        std::cout << P[i] << " ";
    std::cout << "\n";
}

void Solution::printT()
{
    std::cout << "\n=== T (Global): ===\n";
    for (int i = 0; i < amountOfNodes; i++)
        std::cout << T[i] << " ";
    std::cout << "\n";
}

void Solution::solve()
{
    T = gaussElimination(H, P, amountOfNodes);
}

void Element::printH()
{
    std::cout << "\n=== H (Local): ===\n";
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            std::cout << H[i][j] << " ";
        std::cout << "\n";
    }
}

void Element::printHbc()
{
    std::cout << "\n=== Hbc (Local): ===\n";
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            std::cout << Hbc[i][j] << " ";
        std::cout << "\n";
    }
}

void Element::printP()
{
    std::cout << "\n=== P (Local): ===\n";
    for (int i = 0; i < 4; i++)
        std::cout << P[i] << " ";
    std::cout << "\n";
}