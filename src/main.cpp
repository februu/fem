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

  UniversalElement uE;
  uE.initialize();

  Grid grid = globalData.grid;
  Solution solution(grid.amountOfNodes, globalData.initialTemperature, globalData.simulationStep);

  // Loop over all elements
  for (int elementIndex = 0; elementIndex < grid.amountOfElements; elementIndex++)
  {

    Element *element = &(grid.elements[elementIndex]);

    for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
      for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
      {

        // Calculate correct point index
        int currentIndex = i * NUMBER_OF_INTEGRATION_POINTS + j;

        // Calculate jacobian matrix
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

        // Add to Local C
        for (int k = 0; k < 4; k++)
          for (int l = 0; l < 4; l++)
            element->C[k][l] += globalData.density * globalData.specificHeat * uE.N[currentIndex][k] * uE.N[currentIndex][l] * element->jacobians[currentIndex].detJ * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][i] * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][j];
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
              element->Hbc[k][l] += globalData.alfa * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][j] * uE.surfaces[i].N[j][k] * uE.surfaces[i].N[j][l] * detJ;

        // Add to Local P
        for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
          for (int k = 0; k < 4; k++)
            element->P[k] += globalData.alfa * globalData.tot * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][j] * uE.surfaces[i].N[j][k] * detJ;
      }

    // Add local H, Hbc to global H
    // Add local C to global C
    for (int k = 0; k < 4; k++)
      for (int l = 0; l < 4; l++)
      {
        solution.H[element->nodeIds[k]][element->nodeIds[l]] += element->H[k][l];
        solution.H[element->nodeIds[k]][element->nodeIds[l]] += element->Hbc[k][l];
        solution.C[element->nodeIds[k]][element->nodeIds[l]] += element->C[k][l];
      }

    // Add local P to global P
    for (int k = 0; k < 4; k++)
      solution.P[element->nodeIds[k]] += element->P[k];

    // std::cout << "\n========== Element " << elementIndex + 1 << " ========== \n";
    // element->printH();
    // element->printC();
    // element->printHbc();
    // element->printP();
  }

  // std::cout << "\n========== Universal Element ========== \n";
  // uE.print();

  // std::cout << "\n========== Solution ========== \n";
  // solution.printH();
  // solution.printP();
  // solution.printC();

  for (int time = 0; time < globalData.simulationTime; time += globalData.simulationStep)
  {
    // Solve using Gaussian Elimination
    solution.solve();
    std::cout << time + globalData.simulationStep << "\t";
    solution.printTMinMax();
    // solution.printT();
  }
  return 0;
}