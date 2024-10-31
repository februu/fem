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
  Node n1(0.0, 0.0);
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

  for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
  {

    element->jacobians[i].J[0][0] = uE.dN_dKsi[i][0] * grid.nodes[0].x + uE.dN_dKsi[i][1] * grid.nodes[1].x + uE.dN_dKsi[i][2] * grid.nodes[2].x + uE.dN_dKsi[i][3] * grid.nodes[3].x;
    element->jacobians[i].J[0][1] = uE.dN_dKsi[i][0] * grid.nodes[0].y + uE.dN_dKsi[i][1] * grid.nodes[1].y + uE.dN_dKsi[i][2] * grid.nodes[2].y + uE.dN_dKsi[i][3] * grid.nodes[3].y;
    element->jacobians[i].J[1][0] = uE.dN_dEta[i][0] * grid.nodes[0].x + uE.dN_dEta[i][1] * grid.nodes[1].x + uE.dN_dEta[i][2] * grid.nodes[2].x + uE.dN_dEta[i][3] * grid.nodes[3].x;
    element->jacobians[i].J[1][1] = uE.dN_dEta[i][0] * grid.nodes[0].y + uE.dN_dEta[i][1] * grid.nodes[1].y + uE.dN_dEta[i][2] * grid.nodes[2].y + uE.dN_dEta[i][3] * grid.nodes[3].y;
    element->jacobians[i].inverse();

    // Calculate dN_dx and dN_dy
    for (int j = 0; j < 4; j++)
    {
      element->dN_dx[i][j] = element->jacobians[i].invJ[0][0] * uE.dN_dKsi[i][j] + element->jacobians[i].invJ[0][1] * uE.dN_dEta[i][j];
      element->dN_dy[i][j] = element->jacobians[i].invJ[1][0] * uE.dN_dKsi[i][j] + element->jacobians[i].invJ[1][1] * uE.dN_dEta[i][j];
    }

    // FIXME: Calculate H (multiply by weight!!!!)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 4; k++)
        element->H[j][k] += 30 * (element->dN_dx[i][j] * element->dN_dx[i][k] + element->dN_dy[i][j] * element->dN_dy[i][k]) * element->jacobians[i].detJ * 1 * 1;
  }

  // Print H
  for (int j = 0; j < 4; j++)
  {
    for (int k = 0; k < 4; k++)
      std::cout << element->H[j][k] << " ";
    std::cout << "\n";
  }

  return 0;
}