#include "femCalculator.h"

void buildGlobalMatrices(int methodG, GlobalData& gData, Grid& grid,MegaMatrix* megaMatrix_H, MegaMatrix* megaMatrix_C, MegaMatrix* megaMatrix_HBC, std::vector<double>& megaVector) {
    int numElements = (int)gData.data[9];

    ElemUniv elemuniv = elementUniv(methodG);
    ElemUnivC elemunivC = elementUnivC(methodG);

    Surface surfaces[4] = {
        Surface(methodG, 0), Surface(methodG, 1),
        Surface(methodG, 2), Surface(methodG, 3)
    };
    for (int j = 0; j < numElements; j++) {
        for (int i = 0; i < 4; i++) {
            if (grid.elements[j].bc[i] == true) {
                calculateBC(methodG, surfaces[i], i, j, grid, gData.data[3]);
                calculateVectorBC(methodG, surfaces[i], i, j, grid, gData.data[3], gData.data[4]);
            }
        }
    }
    double jtab1[4] = { 0 };
    double jtab2[4] = { 0 };

    for (int i = 0; i < numElements; i++) {
        jtab1[0] = grid.points[grid.elements[i].bl - 1].x; jtab1[1] = grid.points[grid.elements[i].br - 1].x;
        jtab1[2] = grid.points[grid.elements[i].tr - 1].x; jtab1[3] = grid.points[grid.elements[i].tl - 1].x;

        jtab2[0] = grid.points[grid.elements[i].bl - 1].y; jtab2[1] = grid.points[grid.elements[i].br - 1].y;
        jtab2[2] = grid.points[grid.elements[i].tr - 1].y; jtab2[3] = grid.points[grid.elements[i].tl - 1].y;

        Jakobian jakobi = calculateJakobi(methodG, jtab1, jtab2, elemuniv);
        calculate_H_C(gData, elemuniv, elemunivC, jakobi, jtab1, jtab2, gData.data[2], grid, i, surfaces);
    }
}

void runIntegration(const GlobalData& gData, MegaMatrix* megaMatrix_H, MegaMatrix* megaMatrix_C, const std::vector<double>& megaVector) {
    int numNodes = (int)gData.data[8];
    vector<double> t(numNodes, gData.data[5]);

    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            megaMatrix_H->matrix[i][j] += megaMatrix_C->matrix[i][j] / gData.data[1];
        }
    }

    for (double time = 0; time < gData.data[0]; time += gData.data[1]) {
        MegaMatrix matrixtemp(numNodes);
        vector<double> vectorTemp = megaVector;

        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                matrixtemp.matrix[i][j] = megaMatrix_H->matrix[i][j];
            }
        }

        for (int i = 0; i < numNodes; ++i) {
            for (int j = 0; j < numNodes; ++j) {
                vectorTemp[i] += (megaMatrix_C->matrix[i][j] / gData.data[1]) * t[j];
            }
        }

        solveLinearSystem(matrixtemp.matrix, vectorTemp, t, numNodes);

        double min_t = t[0], max_t = t[0];
        for (int i = 0; i < numNodes; i++) {
            if (t[i] > max_t) max_t = t[i];
            if (t[i] < min_t) min_t = t[i];
        }
        cout << "Time: " << time + gData.data[1] << " | Max T: " << max_t << " | Min T: " << min_t << endl;
    }
}