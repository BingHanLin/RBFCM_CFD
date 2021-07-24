#ifndef CONSTANTVALUEIC_HPP
#define CONSTANTVALUEIC_HPP

#include "InitialCondition.hpp"

class ConstantValueIC : public InitialCondition
{
   public:
    ConstantValueIC(const double constValue);

    void fillVector(const size_t nodeID, Eigen::VectorXd& vec) const override;

   private:
    double constValue_;
};

#endif
