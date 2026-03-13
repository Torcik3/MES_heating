#include <iostream>
#include "src/include/mesh.h"
#include "src/include/femCalculator.h"

int main() {
    const string file = "Test1_4_4.txt";
    int methodG = 1;

    vector<double> labels = read_labels(file, 10);
    vector<int> bc_data = read_bc(file);
    GlobalData gData(labels, 4, bc_data);

    int numNodes = (int)gData.data[8];
    int numElements = (int)gData.data[9];

    Grid grid(gData, read_elements(file, numElements), read_points(file, numNodes), bc_data);

    megaMatrix_H = new MegaMatrix(numNodes);
    megaMatrix_HBC = new MegaMatrix(numNodes);
    megaMatrix_C = new MegaMatrix(numNodes);
    megaVector.assign(numNodes, 0.0);

    buildGlobalMatrices(methodG, gData, grid, megaMatrix_H, megaMatrix_C, megaMatrix_HBC, megaVector);
    runIntegration(gData, megaMatrix_H, megaMatrix_C, megaVector);

    delete megaMatrix_H;
    delete megaMatrix_HBC;
    delete megaMatrix_C;

    return 0;
}