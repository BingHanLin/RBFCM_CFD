#ifndef SIMULATIONDOMAIN_HPP
#define SIMULATIONDOMAIN_HPP

#include "json.h"
#include <vector>

/*************************************************************************


*************************************************************************/
template <typename meshType, typename RBFBasisType>
class SimulationDomain
{
   private:
    meshType& myMesh_;
    RBFBasisType& myRBFBasis_;
    nlohmann::json& myParam_;

    std::vector<std::vector<double> > allNodes_;
    int numInn_, numBou_, numAll_;
    double viscous_;
    double delta_t;
    double crankNicolsonEpsilon_;
    int crankNicolsonMaxIter_;

   public:
    SimulationDomain(meshType& mesh, RBFBasisType& RBFBasis,
                     nlohmann::json& param);
    ~SimulationDomain(){};

    void setSimulationType();
};
#endif