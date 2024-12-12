#include "../include/file_parser.h"

void readValue(std::ifstream &file, int &value)
{
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    std::string temp;

    while (iss >> temp)
        if (iss.peek() == std::char_traits<char>::eof())
            value = std::stod(temp);
}

void GlobalData::parseFile(const std::string &fileName)
{
    std::ifstream file("grids/" + fileName);
    if (!file.is_open())
        throw std::runtime_error("Could not open file grids/" + fileName);

    // Loads Global Data
    readValue(file, simulationTime);
    readValue(file, simulationStep);
    readValue(file, conductivity);
    readValue(file, alfa);
    readValue(file, tot);
    readValue(file, initialTemperature);
    readValue(file, density);
    readValue(file, specificHeat);
    readValue(file, grid.amountOfNodes);
    readValue(file, grid.amountOfElements);

    // Loads Nodes
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    grid.nodes = new Node[grid.amountOfNodes];
    for (int i = 0; i < grid.amountOfNodes; i++)
    {
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);
        std::string temp;

        iss.ignore(std::numeric_limits<std::streamsize>::max(), ',');
        iss >> temp;
        grid.nodes[i].x = std::stod(temp);
        iss >> temp;
        grid.nodes[i].y = std::stod(temp);
    }

    // Loads Elements
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    grid.elements = new Element[grid.amountOfElements];
    for (int i = 0; i < grid.amountOfElements; i++)
    {
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);
        std::string temp;

        iss.ignore(std::numeric_limits<std::streamsize>::max(), ',');
        for (int j = 0; j < 4; j++)
        {
            iss >> temp;
            grid.elements[i].nodeIds[j] = std::stoi(temp) - 1;
        }
    }

    // Loads Boundary Conditions
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string value;
        while (iss >> value)
            grid.nodes[std::stoi(value) - 1].isBoundary = true;
    }

    file.close();
}

void GlobalData::print()
{
    std::cout << "\n=== Global Data: ===\n";
    std::cout << "Simulation Time: " << simulationTime << std::endl;
    std::cout << "Simulation Step Time: " << simulationStep << std::endl;
    std::cout << "Conductivity: " << conductivity << std::endl;
    std::cout << "Alfa: " << alfa << std::endl;
    std::cout << "Tot: " << tot << std::endl;
    std::cout << "Initial Temperature: " << initialTemperature << std::endl;
    std::cout << "Density: " << density << std::endl;
    std::cout << "Specific Heat: " << specificHeat << std::endl;
    std::cout << "Number of Nodes: " << grid.amountOfNodes << std::endl;
    std::cout << "Number of Elements: " << grid.amountOfElements << std::endl;

    grid.print();
}