#include <iostream>
#include "../include/file_parser.h"
#include "../include/data_containers.h"
#include "../include/quad.h"

int main()
{

  GlobalData globalData;
  globalData.parseFile("Test2_4_4_MixGrid.txt");
  // globalData.print();

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

    // Add to global H
    for (int k = 0; k < 4; k++)
      for (int l = 0; l < 4; l++)
        solution.H[element->nodeIds[k]][element->nodeIds[l]] += element->H[k][l];
  }

  solution.printH();

  return 0;
}