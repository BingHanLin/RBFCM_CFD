#include "simulationDomain.hpp"
#include "MQBasis2D.hpp"
#include "rectangle.hpp"

template <typename meshType, typename RBFBasisType>
SimulationDomain<meshType, RBFBasisType>::SimulationDomain(
    meshType& mesh, RBFBasisType& RBFBasis, nlohmann::json& param)
    : myMesh_(mesh),
      myRBFBasis_(RBFBasis),
      myParam_(param),
      allNodes_(),
      numInn_(0),
      numBou_(0),
      numAll_(0),
      viscous_(0.0),
      delta_t(0.0),
      crankNicolsonEpsilon_(0.0),
      crankNicolsonMaxIter_(10)
{
    setSimulationType();
}

template <typename meshType, typename RBFBasisType>
void SimulationDomain<meshType, RBFBasisType>::setSimulationType()
{
}

// At end of file
template class SimulationDomain<Rectangle,
                                MQBasis2D>;  // Explicit instantiation