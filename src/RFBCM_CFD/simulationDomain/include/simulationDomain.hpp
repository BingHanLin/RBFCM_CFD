#ifndef SIMULATIONDOMAIN_HPP
#define SIMULATIONDOMAIN_HPP
#include "KDTreeTsaiAdaptor.hpp"
#include "MQBasis.hpp"
#include "json.h"
#include "meshData.hpp"

#include <Eigen/Sparse>
#include <vector>

/*************************************************************************


*************************************************************************/
class SimulationDomain
{
   public:
    SimulationDomain(std::shared_ptr<MeshData> mesh,
                     std::shared_ptr<MQBasis> RBFBasis, nlohmann::json& param);
    ~SimulationDomain(){};

    void showSummary();
    void solveDomain();

    // void exportData();

   private:
    std::shared_ptr<MeshData> myMesh_;
    std::shared_ptr<MQBasis> myRBFBasis_;
    nlohmann::json& myParams_;
    solverTypeEnum solverType_;

    int neighborNum_;
    double endTime_;
    double tStepSize_;
    double density_;
    double viscous_;
    int crankNicolsonMaxIter_;
    double crankNicolsonEpsilon_;

    KDTreeTsaiAdaptor<std::vector<std::vector<double> >, double, 2> kdTree_;

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