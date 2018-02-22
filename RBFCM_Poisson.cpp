#include <iostream>

#include<iostream>

#include <gmm/gmm.h>

#include "GMMRectangular.h"
#include "GMMLrbfcmPoisson.h"
#include "GMMMQBasis2D.h"
using namespace std;
using namespace gmm;

int main()
{

  int m = 20;
  int n = 20;
  double Lx = 1;
  double Ly = 1;
  
  vector< vector<double> > matrix;
  matrix.resize(2, std::vector<double>(3, 0.0));
  matrix[0][1] = 4;
  matrix[1][2] = 5.0;

  vector<double> VX(2),VNX(2);
  VX[0] = 0;
  VX[1] = 1;
  cout << matrix.size()<<endl;
  // Create Mesh
  GMMRECTANGLE MESH(m, n, Lx, Ly);
  // Chose Basis
  GMMMQBasis2D RBFBasis;

  // Set BC condition
  vector<double> Bvalue(2*m+2*n);
  for (int i=0; i<2*m+2*n; i++)
  {
    if (i < m )
      Bvalue[i] = 100;
    else
      Bvalue[i] = 0;   
  }

  // Create Model 
  GMMLrbfcmPoisson<GMMRECTANGLE, GMMMQBasis2D> MODEL(MESH, RBFBasis, Bvalue);
  
  MODEL.SolvePhi();

  return 0;
};