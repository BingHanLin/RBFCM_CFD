
#ifndef GMMRBFCMPOISSON_H
#define GMMRBFCMPOISSON_H

#include<stdio.h>
#include<assert.h>
#include<iostream>
#include<math.h>
#include<vector>
#include <gmm/gmm.h>
#include "KDTreeTsaiAdaptor.h"
#include "Collocation2D.h"
#include "GMMMQBasis2D.h"

using namespace std;
using namespace gmm;
/*************************************************************************


*************************************************************************/
enum BCType_set
     {};

template <typename MeshType, typename RbfBasisType>
class GMMLrbfcmPoisson
{
    private:
        MeshType MESH;
        RbfBasisType RbfBasis;
        vector<vector<double>> AllNodes;
        int Nall, Nint, Nbou;
        vector<double> Phi, Bvalue;
        KDTreeTsaiAdaptor< vector<vector<double> >, double, 2> kdtree;
        
        gmm::row_matrix<rsvector<double>> SystemPhiMatrix; 
        
        void assembly();

    public:
        GMMLrbfcmPoisson (MeshType& MMESH, RbfBasisType& RRbfBasis, const vector<double> BBvalue);
        ~GMMLrbfcmPoisson() {};
        void PrintData();
        void SolvePhi();

};

template <typename MeshType, typename RbfBasisType>
GMMLrbfcmPoisson<MeshType, RbfBasisType>::
GMMLrbfcmPoisson(MeshType& MMESH, RbfBasisType& RRbfBasis, const vector<double> BBvalue)
:MESH(MMESH),RbfBasis(RRbfBasis),Bvalue(BBvalue)
{

    Nall=MESH.Nall;
    Nint=MESH.Nint;
    Nbou=MESH.Nbou;

    AllNodes.resize(Nall);
    vector<double> VX(2);

    for(int i=0; i<Nall; i++)
    {
        VX=MESH.GetAllNode(i+1);
        AllNodes[i].resize(2);
        AllNodes[i][0]=VX[0];
        AllNodes[i][1]=VX[1];
    }

    kdtree = KDTreeTsaiAdaptor<vector<vector<double>>, double, 2>(AllNodes);

    assembly();

    cout<<"Poisson Model is created"<<endl;
}


template <typename MeshType, typename RbfBasisType>
void GMMLrbfcmPoisson<MeshType, RbfBasisType>::
assembly()
{
    int near_num = 5;
    dense_matrix<double> nodes_cloud(near_num,2);
    vector<size_t> neighbours(near_num); // using size_t instead of UINT for compatibility reason with nanoflann.hpp
    vector<double> out_dists_sqr(near_num);
    vector<double> LocalVector;

    gmm::resize(SystemPhiMatrix, Nall, Nall); gmm::clear(SystemPhiMatrix);

    gmm::row_matrix<rsvector<double> > A1(Nint,Nall);
    gmm::resize(A1, Nint, Nall); gmm::clear(A1);
    //==================================
    // go through all interior nodes
    //==================================
    for(int i=0; i<Nint; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdtree.query(i,near_num,&neighbours[0],&out_dists_sqr[0]);
        
        // save nodes cloud in vector
        for(int j=0; j<near_num; j++)
        {
            nodes_cloud(j,0) = AllNodes[neighbours[j]][0];
            nodes_cloud(j,1) = AllNodes[neighbours[j]][1];
        }
    
        gmm::clear(LocalVector);

        LocalVector = Collocation2D(Laplace, nodes_cloud, RbfBasis);
    
        for(int j=0; j<near_num; j++)
        {
            A1(i,neighbours[j]) = LocalVector[j];
        }
    }

    copy(A1, gmm::sub_matrix(SystemPhiMatrix, gmm::sub_interval(0, Nint),
                                              gmm::sub_interval(0, Nall)));


    gmm::resize(A1, Nbou, Nall); gmm::clear(A1);
    //==================================
    // go through all boundary nodes
    //==================================
    for(int i=0; i<Nbou; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdtree.query(i+Nint,near_num,&neighbours[0],&out_dists_sqr[0]);
        
        // save nodes cloud in vector
        for(int j=0; j<near_num; j++)
        {
            nodes_cloud(j,0) = AllNodes[neighbours[j]][0];
            nodes_cloud(j,1) = AllNodes[neighbours[j]][1];
        }
        
        gmm::clear(LocalVector);

        LocalVector = Collocation2D(No_operation, nodes_cloud, RbfBasis);
    
        for(int j=0; j<near_num; j++)
        {
            A1(i,neighbours[j]) = LocalVector[j];
        }
    }

    copy(A1, gmm::sub_matrix(SystemPhiMatrix, gmm::sub_interval(Nint, Nbou),
                                              gmm::sub_interval(0, Nall)));
}


template <typename MeshType, typename RbfBasisType>
void GMMLrbfcmPoisson<MeshType, RbfBasisType>::
SolvePhi()
{
    vector<double> PhiNew(Nall), B1(Nall);
    double rcond;

    Phi.resize(Nall);

    for (int i=0; i<Nall; i++)
    {
        if (i < Nint)
            B1[i] = 0;
        else
            B1[i] = Bvalue[i-Nint];
    }

    lu_solve(SystemPhiMatrix, Phi, B1); //Ax = b
    PrintData();
}

template <typename MeshType, typename RbfBasisType>
void GMMLrbfcmPoisson<MeshType, RbfBasisType>::
PrintData()
{
    char fileoutput[256]="output.txt";

    ofstream fout(fileoutput);

    vector<double> VX(2);

    for(int i=0; i<Nall; i++)
    {
        VX=MESH.GetAllNode(i+1);

        fout<<VX[0]<<" "<<VX[1]<<" "<<Phi[i]<<endl;
    }

    fout.close();

}

#endif