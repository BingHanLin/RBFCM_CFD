#ifndef CONSTANTVALUEBC_HPP
#define CONSTANTVALUEBC_HPP

#include "boundaryCondition.hpp"

class meshData;

class ConstantValueBC : public BoundaryCondition
{
   public:
    ConstantValueBC(const double constValue_);
    // addNode(const int index);

    boundaryConditionType type() override;

   private:
    // std::vector<int> nodeIndices_;
    double constValue_;
};

#endif