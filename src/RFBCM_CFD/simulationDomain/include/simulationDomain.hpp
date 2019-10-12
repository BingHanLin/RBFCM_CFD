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

    double endTime_;
    double tStep_;
    double crankNicolsonEpsilon_;
    int crankNicolsonMaxIter_;
    double density_;
    double viscous_;

    void setUpSimulation();

   public:
    SimulationDomain(meshType& mesh, RBFBasisType& RBFBasis,
                     nlohmann::json& param);
    ~SimulationDomain(){};

    void showSummary();
};
#endif