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
    SimulationDomain(std::shared_ptr<MeshData> mesh,
                     std::shared_ptr<MQBasis> RBFBasis);
    ~SimulationDomain(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    controlData* controlData_;
    std::shared_ptr<MeshData> myMesh_;
    std::shared_ptr<MQBasis> myRBFBasis_;
    solverTypeEnum solverType_;

    int neighborNum_;
    double endTime_;
    double tStepSize_;
    double density_;
    double viscous_;
    int crankNicolsonMaxIter_;
    double crankNicolsonEpsilon_;

    Eigen::SparseMatrix<double> varCoeffMatrix_;
    Eigen::VectorXd varRhs_;
    Eigen::VectorXd varSol_;

    Eigen::SparseMatrix<double> velCoeffMatrix_;
    Eigen::SparseMatrix<double> pressureCoeffMatrix_;

    void setUpSimulation();
    void assembleCoeffMatrix();
    // void assembleRhsMatrix()
};
#endif