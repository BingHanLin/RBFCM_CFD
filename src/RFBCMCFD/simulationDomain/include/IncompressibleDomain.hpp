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
class PUConditionPool;

class IncompressibleDomain
{
   public:
    IncompressibleDomain(ControlData* controlData, MQBasis* RBFBasis,
                         MeshData* meshData, PUConditionPool* conditionPool);
    ~IncompressibleDomain(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    ControlData* controlData_;
    MQBasis* RBFBasis_;
    MeshData* meshData_;
    PUConditionPool* conditionPool_;

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
    Eigen::VectorXd velRhs_;
    Eigen::VectorXd velSol_;
    Eigen::VectorXd preVelSol_;
    Eigen::VectorXd pSol_;
    Eigen::VectorXd prePSol_;

    Eigen::SparseMatrix<double, Eigen::RowMajor> phiCoeffMatrix_;
    Eigen::SparseMatrix<double, Eigen::RowMajor> velCoeffMatrix_;

    Eigen::SparseMatrix<double, Eigen::RowMajor> laplaceMatrix_;
    Eigen::SparseMatrix<double, Eigen::RowMajor> dxMatrix_;
    Eigen::SparseMatrix<double, Eigen::RowMajor> dyMatrix_;
    Eigen::SparseMatrix<double, Eigen::RowMajor> dzMatrix_;

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