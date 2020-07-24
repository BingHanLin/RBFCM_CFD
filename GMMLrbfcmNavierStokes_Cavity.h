
#ifndef GMMRBFCMPOISSON_H
#define GMMRBFCMPOISSON_H

#include <assert.h>
#include <direct.h>
#include <getfem/bgeot_ftool.h>
#include <getfem/getfem_superlu.h>
#include <gmm/gmm.h>
#include <iomanip>  // setprecision
#include <iostream>
#include <math.h>
#include <sstream>  // stringstream
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#include "LRBF/Collocation2D.h"
#include "LRBF/KDTreeTsaiAdaptor.h"

using namespace std;
using namespace gmm;
/*************************************************************************


*************************************************************************/
template <typename MeshType, typename RBFBasisType>
class GMMLrbfcmNavierStokes
{
   private:
    MeshType MESH;
    RBFBasisType RBFBasis;

    vector<vector<double> > AllNodes;
    size_t Nall, Nint, Nbou;
    double viscous, dt;
    double CrankNicolsonEpsilon;
    size_t CrankNicolsonMaxIter;
    vector<double> Bvalue;

    KDTreeTsaiAdaptor<vector<vector<double> >, double, 2> kdtree;

    vector<double> Veln1, Pn1;
    gmm::row_matrix<rsvector<double> > SystemPhiMatrix, SystemVelMatrix;
    gmm::row_matrix<rsvector<double> > IntDxMatrix, IntDyMatrix, BouDxMatrix,
        BouDyMatrix, IntLAPLACEMatrix, BouLAPLACEMatrix;

    void assembly();
    vector<double> CrankNicolsonU(size_t Istep, const vector<double>& Pn,
                                  const vector<double>& Veln);
    vector<double> InsideCrankNicolsonU(size_t Istep, const vector<double>& Pn,
                                        const vector<double>& Veln,
                                        const vector<double>& Veltemp);
    void SolvePU(size_t Istep, vector<double>& Pn1, vector<double>& Veln1,
                 const vector<double>& Pn, const vector<double>& Vels);
    void SaveData(string FolderName, double SaveTime);
    void DeleteOriginData(string FolderName);

   public:
    GMMLrbfcmNavierStokes(MeshType& MMESH, RBFBasisType& RRBFBasis,
                          double vviscous, double ddt,
                          const vector<double> BBvalue);
    ~GMMLrbfcmNavierStokes(){};

    void InitializeField();
};

template <typename MeshType, typename RBFBasisType>
GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::GMMLrbfcmNavierStokes(
    MeshType& MMESH, RBFBasisType& RRBFBasis, double vviscous, double ddt,
    const vector<double> BBvalue)
    : MESH(MMESH),
      RBFBasis(RRBFBasis),
      viscous(vviscous),
      dt(ddt),
      Bvalue(BBvalue)
{
    Nall = MESH.Nall;
    Nint = MESH.Nint;
    Nbou = MESH.Nbou;

    AllNodes.resize(Nall);
    vector<double> VX(2);

    for (size_t i = 0; i < Nall; i++)
    {
        VX = MESH.GetAllNode(i + 1);
        AllNodes[i].resize(2);
        AllNodes[i][0] = VX[0];
        AllNodes[i][1] = VX[1];
    }

    kdtree = KDTreeTsaiAdaptor<vector<vector<double> >, double, 2>(AllNodes);

    assembly();

    cout << "Navier Stokes Model is created" << endl;
}

template <typename MeshType, typename RBFBasisType>
void GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::InitializeField()
{
    Veln1.resize(MESH.Dim * Nall);
    gmm::clear(Veln1);

    Pn1.resize(Nall);
    gmm::clear(Pn1);
}

template <typename MeshType, typename RBFBasisType>
void GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::assembly()
{
    size_t near_num = 9;
    dense_matrix<double> nodes_cloud(near_num, 2);
    vector<size_t> neighbours(
        near_num);  // using size_t instead of UINT for compatibility reason
                    // with nanoflann.hpp
    vector<double> out_dists_sqr(near_num);
    vector<double> LocalVector;

    gmm::resize(SystemPhiMatrix, Nall, Nall);
    gmm::clear(SystemPhiMatrix);
    gmm::resize(SystemVelMatrix, Nall, Nall);
    gmm::clear(SystemVelMatrix);

    gmm::resize(IntDxMatrix, Nint, Nall);
    gmm::clear(IntDxMatrix);
    gmm::resize(IntDyMatrix, Nint, Nall);
    gmm::clear(IntDyMatrix);
    gmm::resize(IntLAPLACEMatrix, Nint, Nall);
    gmm::clear(IntLAPLACEMatrix);

    gmm::resize(BouDxMatrix, Nbou, Nall);
    gmm::clear(BouDxMatrix);
    gmm::resize(BouDyMatrix, Nbou, Nall);
    gmm::clear(BouDyMatrix);
    gmm::resize(BouLAPLACEMatrix, Nbou, Nall);
    gmm::clear(BouLAPLACEMatrix);

    //====================================================================
    // go through all size_terior nodes
    //====================================================================

    for (size_t i = 0; i < Nint; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdtree.query(i, near_num, &neighbours[0], &out_dists_sqr[0]);

        // save nodes cloud in vector
        for (size_t j = 0; j < near_num; j++)
        {
            nodes_cloud(j, 0) = AllNodes[neighbours[j]][0];
            nodes_cloud(j, 1) = AllNodes[neighbours[j]][1];
        }

        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::LAPLACE);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            IntLAPLACEMatrix(i, neighbours[j]) = LocalVector[j];
        }

        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::PARTIAL_D1);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            IntDxMatrix(i, neighbours[j]) = LocalVector[j];
        }

        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::PARTIAL_D2);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            IntDyMatrix(i, neighbours[j]) = LocalVector[j];
        }
    }

    copy(IntLAPLACEMatrix,
         gmm::sub_matrix(SystemPhiMatrix, gmm::sub_interval(0, Nint),
                         gmm::sub_interval(0, Nall)));

    copy(IntLAPLACEMatrix,
         gmm::sub_matrix(SystemVelMatrix, gmm::sub_interval(0, Nint),
                         gmm::sub_interval(0, Nall)));

    for (size_t i = 0; i < Nint; i++)
    {
        SystemVelMatrix(i, i) -= 2.0 / viscous / dt;
    }

    //====================================================================
    // go through all boundary nodes
    //====================================================================
    gmm::row_matrix<rsvector<double> > TempMatrix1;
    gmm::resize(TempMatrix1, Nbou, Nall);
    gmm::clear(TempMatrix1);

    gmm::row_matrix<rsvector<double> > TempMatrix2;
    gmm::resize(TempMatrix2, Nbou, Nall);
    gmm::clear(TempMatrix2);

    vector<double> NormVec;

    for (size_t i = 0; i < Nbou; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdtree.query(i + Nint, near_num, &neighbours[0], &out_dists_sqr[0]);

        // save nodes cloud in vector
        for (size_t j = 0; j < near_num; j++)
        {
            nodes_cloud(j, 0) = AllNodes[neighbours[j]][0];
            nodes_cloud(j, 1) = AllNodes[neighbours[j]][1];
        }

        if (i == 0)  // setup reference Phi value at the first boundary point
        {
            TempMatrix1(i, i + Nint) = 1.0;
        }
        else
        {
            NormVec = MESH.GetOutNormal(i + 1);

            gmm::clear(LocalVector);
            RBFBasis.SetOperatorStatus(
                RBFBasisType::OperatorType::NEUMANN_OPERATOR, NormVec);
            LocalVector = Collocation2D(nodes_cloud, RBFBasis);

            for (size_t j = 0; j < near_num; j++)
            {
                TempMatrix1(i, neighbours[j]) = LocalVector[j];
            }
        }

        // Bounadry SystemVelMatrix
        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::CONSTANT);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            TempMatrix2(i, neighbours[j]) = LocalVector[j];
        }

        // BouDxMatrix
        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::PARTIAL_D1);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            BouDxMatrix(i, neighbours[j]) = LocalVector[j];
        }

        // BouDyMatrix
        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::PARTIAL_D2);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            BouDyMatrix(i, neighbours[j]) = LocalVector[j];
        }

        // BouLAPLACEMatrix
        gmm::clear(LocalVector);
        RBFBasis.SetOperatorStatus(RBFBasisType::OperatorType::LAPLACE);
        LocalVector = Collocation2D(nodes_cloud, RBFBasis);

        for (size_t j = 0; j < near_num; j++)
        {
            BouLAPLACEMatrix(i, neighbours[j]) = LocalVector[j];
        }
    }

    copy(TempMatrix1,
         gmm::sub_matrix(SystemPhiMatrix, gmm::sub_interval(Nint, Nbou),
                         gmm::sub_interval(0, Nall)));

    copy(TempMatrix2,
         gmm::sub_matrix(SystemVelMatrix, gmm::sub_interval(Nint, Nbou),
                         gmm::sub_interval(0, Nall)));

    // cout << SystemVelMatrix << endl;
}

template <typename MeshType, typename RBFBasisType>
vector<double> GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::CrankNicolsonU(
    size_t Istep, const vector<double>& Pn, const vector<double>& Veln)
{
    vector<double> temp1(Veln), temp2(Veln), temp3(Veln);
    double curr_error;
    size_t IterN = 0;

    do
    {
        temp2 = InsideCrankNicolsonU(Istep, Pn, Veln, temp1);

        gmm::add(temp1, gmm::scaled(temp2, -1), temp3);

        curr_error = gmm::vect_norm1(temp3) / gmm::vect_norm1(temp2);

        temp1 = temp2;

        cout << "CrankNicolson error at " << IterN << " = " << fixed
             << setprecision(10) << curr_error << endl;

        IterN++;

        if (IterN > CrankNicolsonMaxIter) break;

    } while (curr_error > CrankNicolsonEpsilon);

    return temp2;
}

template <typename MeshType, typename RBFBasisType>
vector<double>
GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::InsideCrankNicolsonU(
    size_t Istep, const vector<double>& Pn, const vector<double>& Veln,
    const vector<double>& Veltemp)
{
    vector<double> Velout(Veln), VelHalf(Veln), RHS(Veln);

    gmm::clear(Velout);
    gmm::clear(VelHalf);
    gmm::clear(RHS);

    vector<double> RHS_temp1(RHS), RHS_temp2(RHS);

    gmm::add(gmm::scaled(Veln, 0.5), gmm::scaled(Veltemp, 0.5), VelHalf);

    // u
    gmm::mult(IntDxMatrix, gmm::sub_vector(VelHalf, gmm::sub_interval(0, Nall)),
              gmm::sub_vector(RHS_temp1, gmm::sub_interval(0, Nint)));

    gmm::mult(IntDxMatrix,
              gmm::sub_vector(VelHalf, gmm::sub_interval(Nall, Nall)),
              gmm::sub_vector(RHS_temp1, gmm::sub_interval(Nall, Nint)));

    for (size_t i = 0; i < Nint; i++)
    {
        RHS_temp1[i] = RHS_temp1[i] * VelHalf[i];
        RHS_temp1[Nall + i] = RHS_temp1[Nall + i] * VelHalf[i];
    }

    // v
    gmm::mult(IntDyMatrix, gmm::sub_vector(VelHalf, gmm::sub_interval(0, Nall)),
              gmm::sub_vector(RHS_temp2, gmm::sub_interval(0, Nint)));

    gmm::mult(IntDyMatrix,
              gmm::sub_vector(VelHalf, gmm::sub_interval(Nall, Nall)),
              gmm::sub_vector(RHS_temp2, gmm::sub_interval(Nall, Nint)));

    for (size_t i = 0; i < Nint; i++)
    {
        RHS_temp2[i] = RHS_temp2[i] * VelHalf[Nall + i];
        RHS_temp2[Nall + i] = RHS_temp2[Nall + i] * VelHalf[Nall + i];
    }

    gmm::add(RHS_temp1, RHS_temp2, RHS);

    for (size_t i = 0; i < Nint; i++)
    {
        // VelHalf?
        RHS[i] = 2.0 * RHS[i] / viscous - (2.0 / dt / viscous) * Veln[i];
        RHS[Nall + i] = 2.0 * RHS[Nall + i] / viscous -
                        (2.0 / dt / viscous) * Veln[Nall + i];
    }

    // Grad Pn
    gmm::mult_add(IntDxMatrix, gmm::scaled(Pn, 2.0 / viscous),
                  gmm::sub_vector(RHS, gmm::sub_interval(0, Nint)));
    gmm::mult_add(IntDyMatrix, gmm::scaled(Pn, 2.0 / viscous),
                  gmm::sub_vector(RHS, gmm::sub_interval(Nall, Nint)));

    // // LAPLACE Un
    gmm::mult_add(
        IntLAPLACEMatrix,
        gmm::scaled(gmm::sub_vector(Veln, gmm::sub_interval(0, Nall)), -1),
        gmm::sub_vector(RHS, gmm::sub_interval(0, Nint)));
    gmm::mult_add(
        IntLAPLACEMatrix,
        gmm::scaled(gmm::sub_vector(Veln, gmm::sub_interval(Nall, Nall)), -1),
        gmm::sub_vector(RHS, gmm::sub_interval(Nall, Nint)));

    // boundary part ==========================================
    for (size_t i = 0; i < Nbou; i++)
    {
        if (i < 10)
        {
            RHS[Nint + i] = 1.0;
            RHS[Nall + Nint + i] = 0.0;
        }
        else
        {
            RHS[Nint + i] = 0.0;
            RHS[Nall + Nint + i] = 0.0;
        }
    }

    double rcond;

    gmm::SuperLU_solve(SystemVelMatrix,
                       gmm::sub_vector(Velout, gmm::sub_interval(0, Nall)),
                       gmm::sub_vector(RHS, gmm::sub_interval(0, Nall)),
                       rcond);  // Ax = b
    // cout<<"CrankNicolsonU1 : rcond = "<<rcond<<endl;

    gmm::SuperLU_solve(SystemVelMatrix,
                       gmm::sub_vector(Velout, gmm::sub_interval(Nall, Nall)),
                       gmm::sub_vector(RHS, gmm::sub_interval(Nall, Nall)),
                       rcond);  // Ax = b
    // cout<<"CrankNicolsonU2 : rcond = "<<rcond<<endl;

    return Velout;
}

template <typename MeshType, typename RBFBasisType>
void GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::SolvePU(
    size_t Istep, vector<double>& Pn1, vector<double>& Veln1,
    const vector<double>& Pn, const vector<double>& Vels)
{
    gmm::copy(Pn, Pn1);
    gmm::copy(Vels, Veln1);

    vector<double> Phi(Nall), LAPLACEPhi(Nall), GradPhi(2 * Nall), RHS(Nall);
    gmm::clear(Phi);
    gmm::clear(LAPLACEPhi);
    gmm::clear(GradPhi);
    gmm::clear(RHS);

    gmm::mult_add(IntDxMatrix,
                  gmm::scaled(gmm::sub_vector(Vels, gmm::sub_interval(0, Nall)),
                              1.0 / dt),
                  gmm::sub_vector(RHS, gmm::sub_interval(0, Nint)));
    gmm::mult_add(
        IntDyMatrix,
        gmm::scaled(gmm::sub_vector(Vels, gmm::sub_interval(Nall, Nall)),
                    1.0 / dt),
        gmm::sub_vector(RHS, gmm::sub_interval(0, Nint)));

    for (size_t i = 0; i < Nbou; i++)
    {
        if (i == 0)
            RHS[Nint + i] = 0.0;  // reference Phi value
        else
            RHS[Nint + i] = 0.0;
    }

    double rcond;
    gmm::SuperLU_solve(SystemPhiMatrix, Phi, RHS, rcond);
    // cout<<"SolveP (Phi) : rcond = "<< rcond << endl;

    gmm::add(Phi, Pn1);

    gmm::mult(IntLAPLACEMatrix, gmm::scaled(Phi, -viscous * dt / 2.0),
              gmm::sub_vector(LAPLACEPhi, gmm::sub_interval(0, Nint)));
    gmm::mult(BouLAPLACEMatrix, gmm::scaled(Phi, -viscous * dt / 2.0),
              gmm::sub_vector(LAPLACEPhi, gmm::sub_interval(Nint, Nbou)));

    gmm::add(LAPLACEPhi, Pn1);

    double tempP = 0.0 - Pn1[Nint];

    for (size_t i = 0; i < Nall; i++)
        Pn1[i] += tempP;

    gmm::mult_add(IntDxMatrix, gmm::scaled(Phi, -dt),
                  gmm::sub_vector(GradPhi, gmm::sub_interval(0, Nint)));
    // gmm::mult_add(BouDxMatrix, gmm::scaled(Phi,-dt), gmm::sub_vector(GradPhi,
    // gmm::sub_interval(Nint      , Nbou)));

    gmm::mult_add(IntDyMatrix, gmm::scaled(Phi, -dt),
                  gmm::sub_vector(GradPhi, gmm::sub_interval(Nall, Nint)));
    // gmm::mult_add(BouDyMatrix, gmm::scaled(Phi,-dt), gmm::sub_vector(GradPhi,
    // gmm::sub_interval(Nall+Nint , Nbou)));

    gmm::add(GradPhi, Veln1);
}

template <typename MeshType, typename RBFBasisType>
void GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::RunSimulation()
{
    size_t StepN = 500;
    CrankNicolsonEpsilon = 0.00001;
    CrankNicolsonMaxIter = 20;
    vector<double> Veln(Veln1), Pn(Pn1), Vels(Veln1);
    string SaveFolderName = "OutputData";

    DeleteOriginData(SaveFolderName);

    // time loop star
    for (size_t Istep = 1; Istep <= StepN; Istep++)
    {
        Vels = CrankNicolsonU(Istep, Pn, Veln);
        SolvePU(Istep, Pn1, Veln1, Pn, Vels);

        if (Istep % 20 == 0)
        {
            SaveData(SaveFolderName, Istep * dt);
        }

        Veln = Veln1;
        Pn = Pn1;
    }
}

template <typename MeshType, typename RBFBasisType>
void GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::SaveData(string FolderName,
                                                             double SaveTime)
{
    string FilePath, FolderPath;
    stringstream SaveTimeStream;
    SaveTimeStream << fixed << setprecision(5) << SaveTime;
    string SaveTimeStr = SaveTimeStream.str();

    FilePath = ".\\" + FolderName + "\\" + SaveTimeStr + ".dat";
    FolderPath = ".\\" + FolderName;

    cout << "Save data at time = " + SaveTimeStr << endl;

    char charFolderPath[FolderPath.size() + 1];
    strcpy(charFolderPath, FolderPath.c_str());

    mkdir(charFolderPath);

    // Save Data on Nodes
    ofstream fout(FilePath);

    for (size_t i = 0; i < Nall; i++)
    {
        fout << AllNodes[i][0] << ", " << AllNodes[i][1] << ", " << Veln1[i]
             << ", " << Veln1[Nall + i] << ", " << Pn1[i] << endl;
    }

    fout.close();
}

template <typename MeshType, typename RBFBasisType>
void GMMLrbfcmNavierStokes<MeshType, RBFBasisType>::DeleteOriginData(
    string FolderName)
{
    string FolderPath;
    string DeleteString;

    FolderPath = ".\\" + FolderName;
    DeleteString = "del " + FolderPath + " /f /s /q";

    cout << "Delete Original Data in " << FolderName << endl;

    char charDeleteString[DeleteString.size() + 1];
    strcpy(charDeleteString, DeleteString.c_str());

    system(charDeleteString);
}

#endif