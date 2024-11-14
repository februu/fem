#include <iostream>
#include "../include/file_parser.h"
#include "../include/data_containers.h"
#include "../include/quad.h"

double f1(double x)
{
  return 5 * pow(x, 2) + 3 * x + 6;
}

double f2(double x, double y)
{
  return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
}

int main()
{

  // GlobalData globalData;
  // globalData.parseFile("Test1_4_4.txt");
  // globalData.print();

  Grid grid(4, 1);
  Node n1(0.01, -0.01);
  Node n2(0.025, 0.0);
  Node n3(0.025, 0.025);
  Node n4(0.0, 0.025);
  grid.nodes[0] = n1;
  grid.nodes[1] = n2;
  grid.nodes[2] = n3;
  grid.nodes[3] = n4;

  Element e;
  e.nodeIds[0] = 0;
  e.nodeIds[1] = 1;
  e.nodeIds[2] = 2;
  e.nodeIds[3] = 3;
  grid.elements[0] = e;

  UniversalElement uE;
  uE.initialize();

  Element *element = &grid.elements[0];

  // Wersja 2D
  for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
    for (int j = 0; j < NUMBER_OF_INTEGRATION_POINTS; j++)
    {

      // Calculate correct index
      int currentIndex = i * NUMBER_OF_INTEGRATION_POINTS + j;

      element->jacobians[currentIndex].J[0][0] = uE.dN_dKsi[currentIndex][0] * grid.nodes[0].x + uE.dN_dKsi[currentIndex][1] * grid.nodes[1].x + uE.dN_dKsi[currentIndex][2] * grid.nodes[2].x + uE.dN_dKsi[currentIndex][3] * grid.nodes[3].x;
      element->jacobians[currentIndex].J[0][1] = uE.dN_dKsi[currentIndex][0] * grid.nodes[0].y + uE.dN_dKsi[currentIndex][1] * grid.nodes[1].y + uE.dN_dKsi[currentIndex][2] * grid.nodes[2].y + uE.dN_dKsi[currentIndex][3] * grid.nodes[3].y;
      element->jacobians[currentIndex].J[1][0] = uE.dN_dEta[currentIndex][0] * grid.nodes[0].x + uE.dN_dEta[currentIndex][1] * grid.nodes[1].x + uE.dN_dEta[currentIndex][2] * grid.nodes[2].x + uE.dN_dEta[currentIndex][3] * grid.nodes[3].x;
      element->jacobians[currentIndex].J[1][1] = uE.dN_dEta[currentIndex][0] * grid.nodes[0].y + uE.dN_dEta[currentIndex][1] * grid.nodes[1].y + uE.dN_dEta[currentIndex][2] * grid.nodes[2].y + uE.dN_dEta[currentIndex][3] * grid.nodes[3].y;
      element->jacobians[currentIndex].inverse();

      // Calculate dN_dx and dN_dy
      for (int k = 0; k < 4; k++)
      {
        element->dN_dx[currentIndex][k] = element->jacobians[currentIndex].invJ[0][0] * uE.dN_dKsi[currentIndex][k] + element->jacobians[currentIndex].invJ[0][1] * uE.dN_dEta[currentIndex][k];
        element->dN_dy[currentIndex][k] = element->jacobians[currentIndex].invJ[1][0] * uE.dN_dKsi[currentIndex][k] + element->jacobians[currentIndex].invJ[1][1] * uE.dN_dEta[currentIndex][k];
      }

      // FIXME: Calculate H
      for (int k = 0; k < 4; k++)
        for (int l = 0; l < 4; l++)
          element->H[k][l] += 30 * (element->dN_dx[currentIndex][k] * element->dN_dx[currentIndex][l] + element->dN_dy[currentIndex][k] * element->dN_dy[currentIndex][l]) * element->jacobians[currentIndex].detJ * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][i] * gaussWeights[NUMBER_OF_INTEGRATION_POINTS - 1][j];
    }

  // Print H
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
      std::cout << element->H[i][j] << " ";
    std::cout << "\n";
  }

  return 0;
}