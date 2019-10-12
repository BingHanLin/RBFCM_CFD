#include "simulationDomain.hpp"
#include "MQBasis2D.hpp"
#include "rectangle.hpp"
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
      crankNicolsonMaxIter_(10)
{
    setUpSimulation();
    showSummary();
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

    if ("NavierStokes" == myParam_["PhysicsControl"].at("Type"))
    {
        viscous_ = myParam_["PhysicsControl"].at("Viscosity");
        density_ = myParam_["PhysicsControl"].at("Density");
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
    std::cout << std::setfill(' ') << std::setw(56) << "Summary of simulation"
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
}

// explicit instantiation, put this at end of file
template class SimulationDomain<Rectangle, MQBasis2D>;