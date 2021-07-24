#ifndef SIMULATIONDOMAIN_HPP
#define SIMULATIONDOMAIN_HPP
#include "json.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

/*************************************************************************


*************************************************************************/
class MQBasis;
class MeshData;
class ControlData;
class ScalarBCPool;
class ScalarICPool;

class SimulationFlow
{
   public:
    SimulationFlow(ControlData* controlData, MQBasis* RBFBasis,
                   MeshData* meshData, ScalarBCPool* BCPool,
                   ScalarICPool* ICPool);
    ~SimulationFlow(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    ControlData* controlData_;
    MQBasis* RBFBasis_;
    MeshData* meshData_;
    ScalarBCPool* BCPool_;
    ScalarICPool* ICPool_;

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
    void initializeField();
    void setupLinearSystem();
    void assembleCoeffMatrix();
    void assembleRhs();
    void solveMatrix();
    void writeDataToVTK() const;
    void clearVTKDirectory() const;
};
#endif