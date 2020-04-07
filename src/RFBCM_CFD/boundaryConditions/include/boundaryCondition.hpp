#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

#include "enumMap.hpp"

class BoundaryCondition
{
   public:
    BoundaryCondition(){};
    ~BoundaryCondition(){};

    virtual boundaryConditionType type() = 0;

   protected:
};

#endif
