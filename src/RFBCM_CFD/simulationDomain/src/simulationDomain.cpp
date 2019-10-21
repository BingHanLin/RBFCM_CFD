#include "simulationDomain.hpp"
#include "MQBasis2D.hpp"
#include "rectangle.hpp"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

template <typename meshType, typename RBFBasisType>
SimulationDomain<meshType, RBFBasisType>::SimulationDomain(
    meshType& mesh, RBFBasisType& RBFBasis, nlohmann::json& param)
    : myMesh_(mesh),
      myRBFBasis_(RBFBasis),
      myParam_(param),
      viscous_(0.0),
      density_(0.0),
      tStep_(0.0),
      endTime_(0.0),
      crankNicolsonEpsilon_(0.001),
      crankNicolsonMaxIter_(10),
      neighborNum_(5),
      kdTree_()
{
    setUpSimulation();
    showSummary();

    kdTree_ = KDTreeTsaiAdaptor<std::vector<std::vector<double>>, double, 2>(
        myMesh_.getNodes());

    assembleMatrix();
}

template <typename meshType, typename RBFBasisType>
void SimulationDomain<meshType, RBFBasisType>::setUpSimulation()
{
    endTime_ = myParam_["SolverConstrol"].at("TimeStepSize");

    if (endTime_ != 0.0)
    {
        tStep_ = myParam_["SolverConstrol"].at("TimeStepSize");
    }

    try
    {
        neighborNum_ = myParam_["SolverConstrol"].at("NeighborNumber");
    }
    catch (nlohmann::json::out_of_range& e)
    {
        std::cout << "NeighborNumber is not defined, use " << neighborNum_
                  << "." << std::endl;
    }

    if ("NavierStokes" == myParam_["PhysicsControl"].at("Type"))
    {
        viscous_ = myParam_["PhysicsControl"].at("Viscosity");
        density_ = myParam_["PhysicsControl"].at("Density");

        try
        {
            crankNicolsonEpsilon_ =
                myParam_["SolverConstrol"].at("CrankNicolsonEpsilon");
        }
        catch (nlohmann::json::out_of_range& e)
        {
            std::cout << "CrankNicolsonEpsilon is not defined, use "
                      << crankNicolsonEpsilon_ << "." << std::endl;
        }

        try
        {
            crankNicolsonMaxIter_ =
                myParam_["SolverConstrol"].at("CrankNicolsonMaxIter");
        }
        catch (nlohmann::json::out_of_range& e)
        {
            std::cout << "CrankNicolsonMaxIter is not defined, use "
                      << crankNicolsonMaxIter_ << "." << std::endl;
        }
    }
    else if ("Poisson" == myParam_["PhysicsControl"].at("Type"))
    {
        viscous_ = 0.0;
        density_ = 1.0;
    }
}

template <typename meshType, typename RBFBasisType>
void SimulationDomain<meshType, RBFBasisType>::showSummary()
{
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
    std::cout << std::setfill(' ') << std::setw(50) << "Summary of simulation"
              << std::endl;
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::setfill(' ')
              << std::endl;

    std::cout << std::setw(4) << "number of nodes: " << myMesh_.numAllNodes_
              << std::endl;

    if ("NavierStokes" == myParam_["PhysicsControl"].at("Type"))
    {
        std::cout << "density:" << std::setw(4) << density_ << std::endl;
        std::cout << "viscosity:" << std::setw(4) << viscous_ << std::endl;
    }
    else if ("Poisson" == myParam_["PhysicsControl"].at("Type"))
    {
    }

    std::cout << "time step size: " << std::setw(4) << tStep_ << std::endl;
    std::cout << "end time: " << std::setw(4) << endTime_ << std::endl;
    std::cout << "crankNicolson epsilon: " << std::setw(4)
              << crankNicolsonEpsilon_ << std::endl;
    std::cout << "crankNicolson max iter: " << std::setw(4)
              << crankNicolsonMaxIter_ << std::endl;
    std::cout << "Neighbor number: " << std::setw(4) << neighborNum_
              << std::endl;
}

template <typename meshType, typename RBFBasisType>
void SimulationDomain<meshType, RBFBasisType>::assembleMatrix()
{
    // using size_t for compatibility reasonwith nanoflann.hpp
    std::vector<std::vector<double>> nodesCloud(neighborNum_);
    std::vector<size_t> neighbours(neighborNum_);
    std::vector<double> outDistSqr(neighborNum_);
    Eigen::VectorXd localVector;

    // go through all interior nodes
    systemVarMatrix.resize(myMesh_.numAllNodes_, myMesh_.numAllNodes_);
    systemVarMatrix.reserve(
        Eigen::VectorXi::Constant(myMesh_.numAllNodes_, neighborNum_));

    for (int i = 0; i < myMesh_.numInnNodes_; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdTree_.query(i, neighborNum_, &neighbours[0], &outDistSqr[0]);

        // store nodes cloud in vector
        for (int j = 0; j < neighborNum_; j++)
        {
            nodesCloud[j] = myMesh_.getNodes()[neighbours[j]];
        }

        localVector = myRBFBasis_.collectOnNodes(
            nodesCloud, RBFBasisType::operatorType::Laplace);

        for (int j = 0; j < neighborNum_; j++)
        {
            systemVarMatrix.insert(i, neighbours[j]) = localVector(j);
        }
    }

    for (int i = myMesh_.numInnNodes_; i < myMesh_.numAllNodes_; i++)
    {
        // use kdtree find indexes of neighbor nodes
        kdTree_.query(i, neighborNum_, &neighbours[0], &outDistSqr[0]);

        // store nodes cloud in vector
        for (int j = 0; j < neighborNum_; j++)
        {
            nodesCloud[j] = myMesh_.getNodes()[neighbours[j]];
        }

        localVector = myRBFBasis_.collectOnNodes(
            nodesCloud, RBFBasisType::operatorType::IdentityOperation);

        for (int j = 0; j < neighborNum_; j++)
        {
            systemVarMatrix.insert(i, neighbours[j]) = localVector(j);
        }
    }

    systemVarMatrix.makeCompressed();
    std::cout << systemVarMatrix << std::endl;
}

// explicit instantiation, put this at end of file
template class SimulationDomain<Rectangle, MQBasis2D>;