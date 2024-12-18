#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

#include "data_containers.h"

struct GlobalData
{
    int simulationTime;
    int simulationStep;
    int conductivity;
    int alfa;
    int tot;
    int initialTemperature;
    int density;
    int specificHeat;
    Grid grid;

    void parseFile(const std::string &fileName);
    void print();
};
