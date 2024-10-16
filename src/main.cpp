#include "../include/file_parser.h"

int main()
{
  GlobalData globalData;
  globalData.parseFile("Test1_4_4.txt");
  globalData.print();
  return 0;
}