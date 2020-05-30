#ifndef SIMULATIONDOMAIN_HPP
#define SIMULATIONDOMAIN_HPP
#include "KDTreeEigenAdaptor.hpp"
#include "MQBasis.hpp"
#include "controlData.hpp"
#include "json.h"
#include "meshData.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

/*************************************************************************


*************************************************************************/
class SimulationDomain
{
   public:
    SimulationDomain(std::shared_ptr<controlData> inControlData,
                     std::shared_ptr<MeshData> mesh,
                     std::shared_ptr<MQBasis> RBFBasis);
    ~SimulationDomain(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    std::shared_ptr<controlData> controlData_;
    std::shared_ptr<MeshData> myMesh_;
    std::shared_ptr<MQBasis> myRBFBasis_;
    solverType solverType_;
    systemSateType systemSateType_;
    int dim_;

    // physicsControl
    double viscous_;
    double density_;
    double diffusionCoeff_;
    std::array<double, 3> convectionVel_;

    // solverConstrol
    int neighborNum_;
    double endTime_;
    double tStepSize_;
    double writeInterval_;
    double theta_;

    double currentTime_;

    // int crankNicolsonMaxIter_;
    // double crankNicolsonEpsilon_;

    Eigen::SparseMatrix<double> varCoeffMatrix_;
    Eigen::VectorXd varRhs_;
    Eigen::VectorXd varSol_;
    Eigen::VectorXd preVarSol_;

    // Eigen::SparseMatrix<double> velCoeffMatrix_;
    // Eigen::SparseMatrix<double> pressureCoeffMatrix_;
    Eigen::SparseMatrix<double> laplaceMatrix_;
    Eigen::SparseMatrix<double> dxMatrix_;
    Eigen::SparseMatrix<double> dyMatrix_;
    Eigen::SparseMatrix<double> dzMatrix_;

    void setupSimulation();
    void setupInitialCondition();
    void setupLinearSystem();
    void assembleCoeffMatrix();
    void assembleRhs();
    void solveMatrix();
    void writeDataToVTK() const;
    void clearVTKDirectory() const;
};
#endif