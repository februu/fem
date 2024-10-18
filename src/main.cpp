#include "../include/file_parser.h"
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

  std::cout << "Gauss1D: (1): " << gauss1D(f1, 1) << std::endl;
  std::cout << "Gauss1D: (2): " << gauss1D(f1, 2) << std::endl;
  std::cout << "Gauss1D: (3): " << gauss1D(f1, 3) << std::endl;
  std::cout << "Gauss2D: (1): " << gauss2D(f2, 1) << std::endl;
  std::cout << "Gauss2D: (2): " << gauss2D(f2, 2) << std::endl;
  std::cout << "Gauss2D: (3): " << gauss2D(f2, 3) << std::endl;

  std::cout << "Rect1D (5): " << rect1D(f1, 5) << std::endl;
  std::cout << "Rect1D (10): " << rect1D(f1, 10) << std::endl;
  std::cout << "Rect1D (20): " << rect1D(f1, 20) << std::endl;
  std::cout << "Rect1D (50): " << rect1D(f1, 50) << std::endl;

  return 0;
}