#include <iostream>
#include <gmm/gmm.h>

#include "GMMRectangular.h"
#include "GMMLrbfcmNavierStokes_Cavity.h"
#include "LRBF/GMMMQBasis2D.h"

using namespace std;
using namespace gmm;

int main()
{

  int m = 10;
  int n = 10;
  double Lx = 1;
  double Ly = 1;
  double viscous = 1;
  double dt = 0.01;

  // Create Mesh
  GMMRECTANGLE MESH(m, n, Lx, Ly);
  // Chose Basis
  GMMMQBasis2D RBFBasis(6.0);

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
  GMMLrbfcmNavierStokes<GMMRECTANGLE, GMMMQBasis2D> MODEL(MESH, RBFBasis, viscous, dt, Bvalue);

  // Initialize Field
  MODEL.InitializeField();

  MODEL.RunSimulation();

  cout <<"complete"<<endl;
  return 0;
};