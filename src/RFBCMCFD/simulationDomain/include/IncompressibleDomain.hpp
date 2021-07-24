#ifndef INCOMPRESSIBLEDOMAIN_HPP
#define INCOMPRESSIBLEDOMAIN_HPP
#include "enumMap.hpp"
#include "json.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

/*************************************************************************


*************************************************************************/
class MQBasis;
class MeshData;
class ControlData;
class ScalarConditionPool;

class IncompressibleDomain
{
   public:
    IncompressibleDomain(ControlData* controlData, MQBasis* RBFBasis,
                         MeshData* meshData,
                         ScalarConditionPool* conditionPool);
    ~IncompressibleDomain(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    ControlData* controlData_;
    MQBasis* RBFBasis_;
    MeshData* meshData_;
    ScalarConditionPool* conditionPool_;

    // physicsControl
    double viscosity_;
    double density_;

    // solverConstrol
    double neighborRadius_;
    double endTime_;
    double tStepSize_;
    double writeInterval_;
    double crankNicolsonEpsilon_;
    int crankNicolsonMaxIter_;
    size_t dim_;
    double currentTime_;

    // matrix
    Eigen::SparseMatrix<double> varCoeffMatrix_;
    Eigen::VectorXd allRhs_;
    Eigen::VectorXd velSol_;
    Eigen::VectorXd preVelSol_;
    Eigen::VectorXd pSol_;
    Eigen::VectorXd prePSol_;

    Eigen::SparseMatrix<double> phiCoeffMatrix_;
    Eigen::SparseMatrix<double> velCoeffMatrix_;

    Eigen::SparseMatrix<double> laplaceMatrix_;
    Eigen::SparseMatrix<double> dxMatrix_;
    Eigen::SparseMatrix<double> dyMatrix_;
    Eigen::SparseMatrix<double> dzMatrix_;

    void setupSimulation();
    void initializeField();
    void setupLinearSystem();
    void assembleCoeffMatrix();
    void assembleRhs();
    void solveMatrix();
    void writeDataToVTK() const;
    void clearVTKDirectory() const;
};
#endif