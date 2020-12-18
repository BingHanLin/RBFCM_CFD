#ifndef INITIALCONDITION_HPP
#define INITIALCONDITION_HPP

#include "dataStructure.hpp"

#include <memory>

#include <Eigen/Dense>

class InitialCondition
{
   public:
    InitialCondition();

    virtual void fillVector(const size_t nodeID,
                            Eigen::VectorXd& vec) const = 0;
};

#endif
