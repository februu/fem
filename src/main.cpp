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
  Node node1(0.0, 0.0);
  Node node2(0.025, 0.0);
  Node node3(0.025, 0.025);
  Node node4(0.0, 0.025);
  grid.nodes[0] = node1;
  grid.nodes[1] = node2;
  grid.nodes[2] = node3;
  grid.nodes[3] = node4;

  Element element;
  element.nodeIds[0] = 0;
  element.nodeIds[1] = 1;
  element.nodeIds[2] = 2;
  element.nodeIds[3] = 3;
  grid.elements[0] = element;

  UniversalElement uE;
  uE.initialize();

  for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
  {

    grid.elements[0].jacobians[i].J[0][0] = uE.dN_dKsi[i][0] * grid.nodes[0].x + uE.dN_dKsi[i][1] * grid.nodes[1].x + uE.dN_dKsi[i][2] * grid.nodes[2].x + uE.dN_dKsi[i][3] * grid.nodes[3].x;
    grid.elements[0].jacobians[i].J[0][1] = uE.dN_dKsi[i][0] * grid.nodes[0].y + uE.dN_dKsi[i][1] * grid.nodes[1].y + uE.dN_dKsi[i][2] * grid.nodes[2].y + uE.dN_dKsi[i][3] * grid.nodes[3].y;
    grid.elements[0].jacobians[i].J[1][0] = uE.dN_dEta[i][0] * grid.nodes[0].x + uE.dN_dEta[i][1] * grid.nodes[1].x + uE.dN_dEta[i][2] * grid.nodes[2].x + uE.dN_dEta[i][3] * grid.nodes[3].x;
    grid.elements[0].jacobians[i].J[1][1] = uE.dN_dEta[i][0] * grid.nodes[0].y + uE.dN_dEta[i][1] * grid.nodes[1].y + uE.dN_dEta[i][2] * grid.nodes[2].y + uE.dN_dEta[i][3] * grid.nodes[3].y;
    grid.elements[0].jacobians[i].inverse();
  }

  for (int i = 0; i < NUMBER_OF_INTEGRATION_POINTS; i++)
  {
    std::cout << "=== PC " << i + 1 << " ===\n";
    grid.elements[0].jacobians[i].printJ();
    grid.elements[0].jacobians[i].printInvJ();
    std::cout << "det: " << grid.elements[0].jacobians[0].detJ << "\n\n";
  }
  return 0;
}