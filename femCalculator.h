#pragma once
#include <vector>
#include "mesh.h"

void buildGlobalMatrices(int methodG, GlobalData& gData, Grid& grid, MegaMatrix* megaMatrix_H, MegaMatrix* megaMatrix_C, MegaMatrix* megaMatrix_HBC, std::vector<double>& megaVector);
void runIntegration(const GlobalData& gData, MegaMatrix* megaMatrix_H, MegaMatrix* megaMatrix_C, const std::vector<double>& megaVector);