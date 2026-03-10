#include "mesh.h"

// --- Definicje zmiennych globalnych ---
MegaMatrix* megaMatrix_H = nullptr;
MegaMatrix* megaMatrix_HBC = nullptr;
MegaMatrix* megaMatrix_C = nullptr;
vector<double> megaVector;

// --- Implementacje funkcji ---

void solveLinearSystem(vector<vector<double>> A, vector<double> b, vector<double>& x, int n) {
    for (int i = 0; i < n; ++i) {
        A[i].push_back(b[i]);
    }
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        for (int j = 0; j <= n; ++j) {
            std::swap(A[i][j], A[maxRow][j]);
        }

        if (fabs(A[i][i]) < 1e-9) {
            throw std::runtime_error("Uklad rownan jest osobliwy lub niejednoznaczny.");
        }

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j <= n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        x[i] = A[i][n] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
}

vector<double> read_labels(string file, int size) {
    ifstream File(file);
    vector<double> numbers;
    string tag;
    double number;
    int c = 0;

    if (File.is_open()) {
        while (File >> tag >> number && c < size) {
            numbers.push_back(number);
            c++;
        }
        File.close();
    }
    return numbers;
}

vector<Node> read_points(string file, int size) {
    ifstream File(file);
    string line;
    vector<Node> points;
    bool saving = 0;
    string comma;
    int c = 0;

    while (getline(File, line) && c < size) {
        istringstream iss(line);
        if (line.find("*Node") != string::npos) {
            saving = 1;
        }
        if (saving) {
            Node node;
            if (iss >> node.index >> comma >> node.x >> comma >> node.y) {
                points.push_back(node);
                c++;
            }
        }
    }
    return points;
}

vector<Element> read_elements(string file, int size) {
    ifstream File(file);
    string line;
    vector<Element> elements;
    bool saving = 0;
    string comma;
    int c = 0;

    while (getline(File, line) && c < size) {
        istringstream iss(line);
        if (line.find("*Element") != string::npos) {
            saving = 1;
        }
        if (saving) {
            Element elem;
            if (iss >> elem.index >> comma >> elem.tr >> comma >> elem.tl >> comma >> elem.bl >> comma >> elem.br) {
                elements.push_back(elem);
                c++;
            }
        }
    }
    return elements;
}

vector<int> read_bc(string file) {
    ifstream File(file);
    string line;
    vector<int> bc_numbers;
    bool in_bc_section = false;

    while (getline(File, line)) {
        if (line.find("*BC") != string::npos) {
            in_bc_section = true;
            continue;
        }
        if (in_bc_section && !line.empty()) {
            stringstream ss(line);
            string number;
            while (getline(ss, number, ',')) {
                number.erase(0, number.find_first_not_of(" \t"));
                number.erase(number.find_last_not_of(" \t") + 1);
                bc_numbers.push_back(stoi(number));
            }
        }
    }
    File.close();
    return bc_numbers;
}

ElemUnivC elementUnivC(int methodG) {
    int nop = pow(methodG + 1, 2);
    ElemUnivC result(nop);
    double wsp[10] = { 0 };
    double eta, xi;
    switch (methodG) {
    case 1:
        wsp[0] = -0.5774;
        wsp[1] = 0.5774;
        break;
    case 2:
        wsp[0] = -sqrt(3.0 / 5.0);
        wsp[1] = 0;
        wsp[2] = sqrt(3.0 / 5.0);
        break;
    }
    int counter = 0;
    for (int h = 0; h < sqrt(nop); h++) {
        eta = wsp[h];
        for (int j = 0; j < sqrt(nop); j++) {
            xi = wsp[j];
            result.n_C[counter][0] = 0.25 * (1 - eta) * (1 - xi);
            result.n_C[counter][1] = 0.25 * (1 - eta) * (1 + xi);
            result.n_C[counter][2] = 0.25 * (1 + eta) * (1 + xi);
            result.n_C[counter][3] = 0.25 * (1 + eta) * (1 - xi);
            counter++;
        }
    }
    return result;
}

ElemUniv elementUniv(int methodG) {
    int nop = pow(methodG + 1, 2);
    ElemUniv result(nop);
    double wsp[3];
    double w[3];

    switch (methodG) {
    case 0:
        wsp[0] = 0; w[0] = 2; break;
    case 1:
        wsp[0] = -1.0 / sqrt(3.0); wsp[1] = 1.0 / sqrt(3.0);
        w[0] = 1.0; w[1] = 1.0; break;
    case 2:
        wsp[0] = -sqrt(3.0 / 5.0); wsp[1] = 0; wsp[2] = sqrt(3.0 / 5.0);
        w[0] = 5.0 / 9.0; w[1] = 8.0 / 9.0; w[2] = 5.0 / 9.0; break;
    }
    int counter = 0;

    for (int h = 0; h < sqrt(nop); h++) {
        double eta = wsp[h];
        for (int j = 0; j < sqrt(nop); j++) {
            double xi = wsp[j];

            result.dN_dxsi[counter][0] = -0.25 * (1 - eta);
            result.dN_dxsi[counter][1] = 0.25 * (1 - eta);
            result.dN_dxsi[counter][2] = 0.25 * (1 + eta);
            result.dN_dxsi[counter][3] = -0.25 * (1 + eta);

            result.dN_deta[counter][0] = -0.25 * (1 - xi);
            result.dN_deta[counter][1] = -0.25 * (1 + xi);
            result.dN_deta[counter][2] = 0.25 * (1 + xi);
            result.dN_deta[counter][3] = 0.25 * (1 - xi);
            counter++;
        }
    }
    return result;
}

Jakobian calculateJakobi(int methodG, double x[], double y[], ElemUniv elemuUniv) {
    int nop = pow(methodG + 1, 2);
    Jakobian result(nop);

    for (int i = 0; i < nop; i++) {
        result.j[i][0] = x[0] * elemuUniv.dN_dxsi[i][0] + x[1] * elemuUniv.dN_dxsi[i][1] + x[2] * elemuUniv.dN_dxsi[i][2] + x[3] * elemuUniv.dN_dxsi[i][3];
        result.j[i][1] = y[0] * elemuUniv.dN_dxsi[i][0] + y[1] * elemuUniv.dN_dxsi[i][1] + y[2] * elemuUniv.dN_dxsi[i][2] + y[3] * elemuUniv.dN_dxsi[i][3];
        result.j[i][2] = x[0] * elemuUniv.dN_deta[i][0] + x[1] * elemuUniv.dN_deta[i][1] + x[2] * elemuUniv.dN_deta[i][2] + x[3] * elemuUniv.dN_deta[i][3];
        result.j[i][3] = y[0] * elemuUniv.dN_deta[i][0] + y[1] * elemuUniv.dN_deta[i][1] + y[2] * elemuUniv.dN_deta[i][2] + y[3] * elemuUniv.dN_deta[i][3];
    }

    for (int j = 0; j < nop; j++) {
        result.detJ[j] = result.j[j][0] * result.j[j][3] - result.j[j][1] * result.j[j][2];
        result.j1[j][0] = result.j[j][3] / result.detJ[j];
        result.j1[j][1] = -result.j[j][1] / result.detJ[j];
        result.j1[j][2] = -result.j[j][2] / result.detJ[j];
        result.j1[j][3] = result.j[j][0] / result.detJ[j];
    }
    return result;
}

Matrix calculateBC(int methodG, Surface surface, int side, int elementIndex, Grid& grid, double alfa) {
    double w[3];
    double h1[4][4] = { 0 };
    double length;
    Matrix result;
    switch (methodG) {
    case 0: w[0] = 2; break;
    case 1: w[0] = 1.0; w[1] = 1.0; break;
    case 2: w[0] = 5.0 / 9.0; w[1] = 8.0 / 9.0; w[2] = 5.0 / 9.0; break;
    }

    switch (side) {
    case 0: length = sqrt(pow(grid.points[grid.elements[elementIndex].br - 1].x - grid.points[grid.elements[elementIndex].bl - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].br - 1].y - grid.points[grid.elements[elementIndex].bl - 1].y, 2)); break;
    case 1: length = sqrt(pow(grid.points[grid.elements[elementIndex].tr - 1].x - grid.points[grid.elements[elementIndex].br - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].tr - 1].y - grid.points[grid.elements[elementIndex].br - 1].y, 2)); break;
    case 2: length = sqrt(pow(grid.points[grid.elements[elementIndex].tr - 1].x - grid.points[grid.elements[elementIndex].tl - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].tr - 1].y - grid.points[grid.elements[elementIndex].tl - 1].y, 2)); break;
    case 3: length = sqrt(pow(grid.points[grid.elements[elementIndex].tl - 1].x - grid.points[grid.elements[elementIndex].bl - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].tl - 1].y - grid.points[grid.elements[elementIndex].bl - 1].y, 2)); break;
    }

    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            h1[j][k] += w[0] * surface.N[0][j] * surface.N[0][k] * alfa;
            h1[j][k] += w[1] * surface.N[1][j] * surface.N[1][k] * alfa;
            h1[j][k] *= length / 2;
            result.matrix[j][k] += h1[j][k];
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            grid.elements[elementIndex].BC[i][j] += h1[i][j];
        }
    }
    return result;
}

void calculateVectorBC(int methodG, Surface surface, int side, int elementIndex, Grid& grid, double alfa, double Tot) {
    double p[4] = { 0 };
    double w[10];
    double length = 0;
    switch (methodG) {
    case 0: w[0] = 2; break;
    case 1: w[0] = 1.0; w[1] = 1.0; break;
    case 2: w[0] = 5.0 / 9.0; w[1] = 8.0 / 9.0; w[2] = 5.0 / 9.0; break;
    }

    switch (side) {
    case 0: length = sqrt(pow(grid.points[grid.elements[elementIndex].bl - 1].x - grid.points[grid.elements[elementIndex].br - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].bl - 1].y - grid.points[grid.elements[elementIndex].br - 1].y, 2)); break;
    case 1: length = sqrt(pow(grid.points[grid.elements[elementIndex].br - 1].x - grid.points[grid.elements[elementIndex].tr - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].br - 1].y - grid.points[grid.elements[elementIndex].tr - 1].y, 2)); break;
    case 2: length = sqrt(pow(grid.points[grid.elements[elementIndex].tr - 1].x - grid.points[grid.elements[elementIndex].tl - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].tr - 1].y - grid.points[grid.elements[elementIndex].tl - 1].y, 2)); break;
    case 3: length = sqrt(pow(grid.points[grid.elements[elementIndex].tl - 1].x - grid.points[grid.elements[elementIndex].bl - 1].x, 2) + pow(grid.points[grid.elements[elementIndex].tl - 1].y - grid.points[grid.elements[elementIndex].bl - 1].y, 2)); break;
    }

    for (int j = 0; j < 4; j++) {
        p[j] += Tot * surface.N[0][j] * w[0];
        p[j] += Tot * surface.N[1][j] * w[1];
        p[j] *= alfa * length / 2;
        grid.elements[elementIndex].p[j] += p[j];
    }
}

void calculate_H_C(GlobalData& GData, ElemUniv elemUniv, ElemUnivC elemUnivC, Jakobian jakobi, double x[], double y[], double wsp_ciep, Grid& grid, int slon, Surface surfaces[]) {
    vector<vector<double>> ndx(GData.npc, vector<double>(4, 0.0));
    vector<vector<double>> ndy(GData.npc, vector<double>(4, 0.0));
    vector<Matrix> matrixH(GData.npc);
    vector<Matrix> matrixC(GData.npc);
    Matrix result;
    Matrix resultC;

    int counter = 0;
    double wsp[10];
    double w[10];
    vector<double> vec(GData.npc);
    int tab[4] = { grid.elements[slon].bl, grid.elements[slon].br, grid.elements[slon].tr, grid.elements[slon].tl };

    switch (GData.npc) {
    case 0: wsp[0] = 0; w[0] = 2; break;
    case 4: wsp[0] = -1.0 / sqrt(3.0); wsp[1] = 1.0 / sqrt(3.0); w[0] = 1.0; w[1] = 1.0; break;
    case 9: wsp[0] = -sqrt(3.0 / 5.0); wsp[1] = 0; wsp[2] = sqrt(3.0 / 5.0); w[0] = 5.0 / 9.0; w[1] = 8.0 / 9.0; w[2] = 5.0 / 9.0; break;
    }

    for (int i = 0; i < sqrt(GData.npc); i++) {
        for (int j = 0; j < sqrt(GData.npc); j++) {
            vec[counter] = w[i] * w[j];
            counter++;
        }
    }

    for (int i = 0; i < GData.npc; i++) {
        for (int h = 0; h < 4; h++) {
            ndx[i][h] = jakobi.j1[i][0] * elemUniv.dN_dxsi[i][h] + jakobi.j1[i][1] * elemUniv.dN_deta[i][h];
            ndy[i][h] = jakobi.j1[i][2] * elemUniv.dN_dxsi[i][h] + jakobi.j1[i][3] * elemUniv.dN_deta[i][h];
        }
    }

    for (int i = 0; i < GData.npc; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                matrixH[i].matrix[j][k] = wsp_ciep * (ndx[i][j] * ndx[i][k] + ndy[i][j] * ndy[i][k]) * jakobi.detJ[i];
                matrixC[i].matrix[j][k] = GData.data[6] * GData.data[7] * jakobi.detJ[i] * elemUnivC.n_C[i][j] * elemUnivC.n_C[i][k];
            }
        }
    }

    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < GData.npc; i++) {
                result.matrix[j][k] += matrixH[i].matrix[j][k] * vec[i];
                resultC.matrix[j][k] += matrixC[i].matrix[j][k] * vec[i];
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.matrix[i][j] += grid.elements[slon].BC[i][j];
        }
    }

    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            megaMatrix_H->matrix[tab[j] - 1][tab[k] - 1] += result.matrix[j][k];
            megaMatrix_C->matrix[tab[j] - 1][tab[k] - 1] += resultC.matrix[j][k];
        }
    }

    for (int j = 0; j < 4; j++) {
        megaVector[tab[j] - 1] += grid.elements[slon].p[j];
    }
}