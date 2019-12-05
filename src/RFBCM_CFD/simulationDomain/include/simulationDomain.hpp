#ifndef SIMULATIONDOMAIN_HPP
#define SIMULATIONDOMAIN_HPP

#include "KDTreeTsaiAdaptor.hpp"
#include "json.h"

#include <Eigen/Sparse>
#include <vector>

/*************************************************************************


*************************************************************************/
template <typename meshType, typename RBFBasisType>
class SimulationDomain
{
   private:
    meshType& myMesh_;
    RBFBasisType& myRBFBasis_;
    nlohmann::json& myParams_;

    int neighborNum_;
    double endTime_;
    double tStep_;
    double crankNicolsonEpsilon_;
    int crankNicolsonMaxIter_;
    double density_;
    double viscous_;

    KDTreeTsaiAdaptor<std::vector<std::vector<double> >, double, 2> kdTree_;

    Eigen::SparseMatrix<double> systemVarMatrix_;
    Eigen::SparseMatrix<double> systemPressureMatrix_;
    Eigen::SparseMatrix<double> systemVelxMatrix_;
    Eigen::SparseMatrix<double> systemVelyMatrix_;
    Eigen::SparseMatrix<double> systemVelzMatrix_;

    Eigen::VectorXd rhs_;
    Eigen::VectorXd solution_;

    void setUpSimulation();
    void assembleMatrix();

   public:
    SimulationDomain(meshType& mesh, RBFBasisType& RBFBasis,
                     nlohmann::json& param);
    ~SimulationDomain(){};

    void exportData();
    void solveDomain();

    void showSummary();
};
#endif