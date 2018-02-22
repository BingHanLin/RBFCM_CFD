using namespace std;
using namespace gmm;
#include "GMMMQBasis2D.h"
template<typename T>
vector<T> Collocation2D(const int Operator2D, const dense_matrix<T> & nodes_cloud, GMMMQBasis2D & RbfBasis)
{
    int DIMENSION = 2;
    // check the input vector
    assert(mat_ncols(nodes_cloud)==DIMENSION);
    
    int Nearnum = mat_nrows(nodes_cloud);

    dense_matrix<T> Phi(Nearnum,Nearnum);
    vector<T> temp(Nearnum), LPhi(Nearnum);
    vector<T> VX1(DIMENSION), VX2(DIMENSION);

    for(int i=0; i<Nearnum; i++)
    {
        for(int j=0; j<Nearnum; j++)
        {
            VX1[0] = nodes_cloud(i,0);
            VX1[1] = nodes_cloud(i,1);

            VX2[0] = nodes_cloud(j,0);
            VX2[1] = nodes_cloud(j,1);

            Phi(j,i) /* transposed */ = RbfBasis.SetOperator(0, VX1, VX2);

            if (i == 0)
                LPhi[j] = RbfBasis.SetOperator(Operator2D, VX1, VX2);
        }
    }

    // ! 考慮之後加上判斷 ill condition 
    lu_solve(Phi, temp, LPhi);
    return temp;
}