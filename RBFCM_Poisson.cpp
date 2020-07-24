#include <gmm/gmm.h>
#include <iostream>

#include "GMMLrbfcmPoisson.h"
#include "GMMRectangular.h"
#include "LRBF/GMMMQBasis2D.h"

using namespace std;
using namespace gmm;

int main()
{
    size_t m = 10;
    size_t n = 10;
    double Lx = 1;
    double Ly = 1;

    // Create Mesh
    GMMRECTANGLE MESH(m, n, Lx, Ly);
    // Chose Basis
    GMMMQBasis2D RBFBasis(1.0);

    // Set BC condition
    vector<double> Bvalue(2 * m + 2 * n);
    for (size_t i = 0; i < 2 * m + 2 * n; i++)
    {
        if (i < m)
            Bvalue[i] = 100;
        else
            Bvalue[i] = 0;
    }

    // Create Model
    GMMLrbfcmPoisson<GMMRECTANGLE, GMMMQBasis2D> MODEL(MESH, RBFBasis, Bvalue);

    MODEL.SolvePhi();

    return 0;
};