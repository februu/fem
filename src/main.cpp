#include <iostream>
#include <iomanip>
#include "../include/file_parser.h"
#include "../include/data_containers.h"
#include "../include/quad.h"

int main()
{

  GlobalData globalData;
  globalData.parseFile("Test1_4_4.txt");
  // globalData.print();

  std::cout << std::setprecision(5);

  UniversalElement uE;
  uE.initialize();

  Grid grid = globalData.grid;
  Solution solution(grid.amountOfNodes);

  // Loop over all elements
  for (int elementIndex = 0; elementIndex < grid.amountOfElements; elementIndex++)
  {

    Element *element = &(grid.elements[elementIndex]);

    for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
      for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
      {

        // Calculate correct index
        int currentIndex = i * NUMBER_OF_INTEGRATION_POINTS + j;

        element->jacobians[currentIndex].J[0][0] = uE.dN_dKsi[currentIndex][0] * grid.nodes[element->nodeIds[0]].x + uE.dN_dKsi[currentIndex][1] * grid.nodes[element->nodeIds[1]].x + uE.dN_dKsi[currentIndex][2] * grid.nodes[element->nodeIds[2]].x + uE.dN_dKsi[currentIndex][3] * grid.nodes[element->nodeIds[3]].x;
        element->jacobians[currentIndex].J[0][1] = uE.dN_dKsi[currentIndex][0] * grid.nodes[element->nodeIds[0]].y + uE.dN_dKsi[currentIndex][1] * grid.nodes[element->nodeIds[1]].y + uE.dN_dKsi[currentIndex][2] * grid.nodes[element->nodeIds[2]].y + uE.dN_dKsi[currentIndex][3] * grid.nodes[element->nodeIds[3]].y;
        element->jacobians[currentIndex].J[1][0] = uE.dN_dEta[currentIndex][0] * grid.nodes[element->nodeIds[0]].x + uE.dN_dEta[currentIndex][1] * grid.nodes[element->nodeIds[1]].x + uE.dN_dEta[currentIndex][2] * grid.nodes[element->nodeIds[2]].x + uE.dN_dEta[currentIndex][3] * grid.nodes[element->nodeIds[3]].x;
        element->jacobians[currentIndex].J[1][1] = uE.dN_dEta[currentIndex][0] * grid.nodes[element->nodeIds[0]].y + uE.dN_dEta[currentIndex][1] * grid.nodes[element->nodeIds[1]].y + uE.dN_dEta[currentIndex][2] * grid.nodes[element->nodeIds[2]].y + uE.dN_dEta[currentIndex][3] * grid.nodes[element->nodeIds[3]].y;
        element->jacobians[currentIndex].inverse();

        // Calculate dN_dx and dN_dy
        for (int k = 0; k < 4; k++)
        {
          element->dN_dx[currentIndex][k] = element->jacobians[currentIndex].invJ[0][0] * uE.dN_dKsi[currentIndex][k] + element->jacobians[currentIndex].invJ[0][1] * uE.dN_dEta[currentIndex][k];
          element->dN_dy[currentIndex][k] = element->jacobians[currentIndex].invJ[1][0] * uE.dN_dKsi[currentIndex][k] + element->jacobians[currentIndex].invJ[1][1] * uE.dN_dEta[currentIndex][k];
        }

        // Add to Local H
        for (int k = 0; k < 4; k++)
          for (int l = 0; l < 4; l++)
            element->H[k][l] += globalData.conductivity * (element->dN_dx[currentIndex][k] * element->dN_dx[currentIndex][l] + element->dN_dy[currentIndex][k] * element->dN_dy[currentIndex][l]) * element->jacobians[currentIndex].detJ * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][i] * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][j];
      }

    // Check each of the 4 sides for BC
    for (int i = 0; i < 4; i++)
      if (grid.nodes[element->nodeIds[i % 4]].isBoundary && grid.nodes[element->nodeIds[(i + 1) % 4]].isBoundary)
      {
        int firstNode = i % 4;
        int secondNode = (i + 1) % 4;

        double sideLength = sqrt(pow(grid.nodes[element->nodeIds[firstNode]].x - grid.nodes[element->nodeIds[secondNode]].x, 2) + pow(grid.nodes[element->nodeIds[firstNode]].y - grid.nodes[element->nodeIds[secondNode]].y, 2));
        double detJ = sideLength / 2;

        // Add to Local Hbc
        for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
          for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
              if (k == firstNode || k == secondNode || l == firstNode || l == secondNode)
                element->Hbc[k][l] += globalData.alfa * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][j] * uE.surfaces[i].N[j][k] * uE.surfaces[i].N[j][l] * detJ;
      }

    // Add to global H
    for (int k = 0; k < 4; k++)
      for (int l = 0; l < 4; l++)
        solution.H[element->nodeIds[k]][element->nodeIds[l]] += element->H[k][l];

    std::cout << "\nElement " << elementIndex + 1 << ":";
    element->printHbc();
  }

  // solution.printH();

  return 0;
}