#ifndef MES_H
#define MES_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std; // Uwaga: w dużych projektach unika się tego w plikach .h, ale zostawiamy dla zgodności z Twoim kodem.

class Node {
public:
    int index;
    double x, y;
};

struct Jakobian {
public:
    vector<vector<double>> j;
    vector<vector<double>> j1;
    vector<double> detJ;
    int npc;

    Jakobian(int npc) : npc(npc),
    j(npc, vector<double>(4, 0.0)),
    j1(npc, vector<double>(4, 0.0)),
    detJ(npc, 0.0) {}
};

class Element {
public:
    int index;
    int tr, tl, bl, br;
    bool bc[4] = { false };
    double BC[4][4] = { 0 };
    double p[4] = { 0 };

    void setBoundaryConditions(const vector<int>& edge) {
        bool test[4] = { false };
        for (size_t i = 0; i < edge.size(); i++) {
            if (this->bl == edge[i]) test[0] = true;
            if (this->br == edge[i]) test[1] = true;
            if (this->tr == edge[i]) test[2] = true;
            if (this->tl == edge[i]) test[3] = true;
        }

        if (test[0] && test[1]) this->bc[0] = true;
        if (test[1] && test[2]) this->bc[1] = true;
        if (test[2] && test[3]) this->bc[2] = true;
        if (test[3] && test[0]) this->bc[3] = true;
    }
};

class GlobalData {
public:
    vector<double> data;
    int npc;
    vector<int> bc;

    GlobalData(const vector<double>& a, int npc1, const vector<int>& bc1) {
        data = a;
        bc = bc1;
        npc = npc1;
    }
};

class Grid {
public:
    int nN;
    int nE;
    vector<Element> elements;
    vector<Node> points;

    Grid(GlobalData gData, const vector<Element>& arrayE, const vector<Node>& arrayP, const vector<int>& edge) {
        nN = gData.data[8];
        nE = gData.data[9];
        elements = arrayE;
        points = arrayP;

        for (int i = 0; i < nE; i++) {
            elements[i].setBoundaryConditions(edge);
        }
    }
};

struct ElemUniv {
public:
    vector<vector<double>> dN_dxsi;
    vector<vector<double>> dN_deta;

    ElemUniv(int npc) : dN_dxsi(npc, vector<double>(4, 0.0)),
                        dN_deta(npc, vector<double>(4, 0.0)) {}
};

struct ElemUnivC {
public:
    vector<vector<double>> n_C;

    ElemUnivC(int npc) : n_C(npc, vector<double>(4, 0.0)) {}
};

struct Matrix {
    double matrix[4][4] = { 0 };
};

struct MegaMatrix {
    vector<vector<double>> matrix;

    MegaMatrix(int size) : matrix(size, vector<double>(size, 0.0)) {}
};

struct Surface {
    vector<vector<double>> N;
    int methodG;

    Surface(int methodG, int side) {
        int nop = pow(methodG + 1, 2);
        this->methodG = methodG;
        N.assign(methodG + 1, vector<double>(4, 0.0));

        double ksi[10][10] = { 0 };
        double eta[10][10] = { 0 };

        if (this->methodG == 1) {
            ksi[0][0] = -1.0 / sqrt(3.0); ksi[0][1] = 1.0 / sqrt(3.0);
            ksi[1][0] = 1; ksi[1][1] = 1;
            ksi[2][0] = 1.0 / sqrt(3.0); ksi[2][1] = -1.0 / sqrt(3.0);
            ksi[3][0] = -1.0; ksi[3][1] = -1.0;

            eta[0][0] = -1.0; eta[0][1] = -1.0;
            eta[1][0] = -1.0 / sqrt(3.0); eta[1][1] = 1.0 / sqrt(3.0);
            eta[2][0] = 1.0; eta[2][1] = 1.0;
            eta[3][0] = 1.0 / sqrt(3.0); eta[3][1] = -1.0 / sqrt(3.0);
        }
        if (this->methodG == 2) {
            ksi[0][0] = -sqrt(3.0 / 5.0); ksi[0][1] = 0; ksi[0][2] = sqrt(3.0 / 5.0);
            ksi[1][0] = 1; ksi[1][1] = 1; ksi[1][2] = 1;
            ksi[2][0] = sqrt(3.0 / 5.0); ksi[2][1] = 0; ksi[2][2] = -sqrt(3.0 / 5.0);
            ksi[3][0] = -1.0; ksi[3][1] = -1.0; ksi[3][2] = -1.0;

            eta[0][0] = -1.0; eta[0][1] = -1.0; eta[0][2] = -1.0;
            eta[1][0] = -sqrt(3.0 / 5.0); eta[1][1] = 0; eta[1][2] = sqrt(3.0 / 5.0);
            eta[2][0] = 1.0; eta[2][1] = 1.0; eta[2][2] = 1.0;
            eta[3][0] = sqrt(3.0 / 5.0); eta[3][1] = 0; eta[3][2] = -sqrt(3.0 / 5.0);
        }

        for (int i = 0; i < methodG + 1; i++) {
            N[i][0] = 0.25 * (1 - ksi[side][i]) * (1 - eta[side][i]);
            N[i][1] = 0.25 * (1 + ksi[side][i]) * (1 - eta[side][i]);
            N[i][2] = 0.25 * (1 + ksi[side][i]) * (1 + eta[side][i]);
            N[i][3] = 0.25 * (1 - ksi[side][i]) * (1 + eta[side][i]);
        }
    }
};

// --- Zmienne globalne (tylko deklaracje dzięki 'extern') ---
extern MegaMatrix* megaMatrix_H;
extern MegaMatrix* megaMatrix_HBC;
extern MegaMatrix* megaMatrix_C;
extern vector<double> megaVector;

// --- Prototypy funkcji ---
void solveLinearSystem(vector<vector<double>> A, vector<double> b, vector<double>& x, int n);
vector<double> read_labels(string file, int size);
vector<Node> read_points(string file, int size);
vector<Element> read_elements(string file, int size);
vector<int> read_bc(string file);
ElemUnivC elementUnivC(int methodG);
ElemUniv elementUniv(int methodG);
Jakobian calculateJakobi(int methodG, double x[], double y[], ElemUniv elemuUniv);
Matrix calculateBC(int methodG, Surface surface, int side, int elementIndex, Grid& grid, double alfa);
void calculateVectorBC(int methodG, Surface surface, int side, int elementIndex, Grid& grid, double alfa, double Tot);
void calculate_H_C(GlobalData& GData, ElemUniv elemUniv, ElemUnivC elemUnivC, Jakobian jakobi, double x[], double y[], double wsp_ciep, Grid& grid, int slon, Surface surfaces[]);

#endif // MES_H