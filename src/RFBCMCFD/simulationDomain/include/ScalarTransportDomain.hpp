#ifndef SIMULATIONDOMAIN_HPP
#define SIMULATIONDOMAIN_HPP
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

class ScalarTransportDomain
{
   public:
    ScalarTransportDomain(ControlData* controlData, MQBasis* RBFBasis,
                          MeshData* meshData,
                          ScalarConditionPool* conditionPool);
    ~ScalarTransportDomain(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    ControlData* controlData_;
    MQBasis* RBFBasis_;
    MeshData* meshData_;
    ScalarConditionPool* conditionPool_;

    // physicsControl
    double density_;
    double diffusionCoeff_;
    std::array<double, 3> convectionVel_;

    // solverConstrol
    double neighborRadius_;
    double endTime_;
    double tStepSize_;
    double writeInterval_;
    double theta_;
    systemSateType systemSateType_;
    size_t dim_;
    double currentTime_;

    // matrix
    Eigen::SparseMatrix<double, Eigen::RowMajor> varCoeffMatrix_;
    Eigen::VectorXd varRhs_;
    Eigen::VectorXd varSol_;
    Eigen::VectorXd preVarSol_;

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