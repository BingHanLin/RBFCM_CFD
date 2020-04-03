#ifndef MQBasis_HPP
#define MQBasis_HPP

#include "dataStructure.hpp"

#include <Eigen/Dense>
#include <vector>

/*************************************************************************


*************************************************************************/
class MQBasis
{
   public:
    enum operatorType
    {
        IdentityOperation,
        Laplace,
        Partial_D1,
        Partial_D2,
        Partial_D3,
        NEUMANN_OPERATOR
    };

    MQBasis(double CCC = 1.0) : cc_(CCC){};

    ~MQBasis(){};

    Eigen::VectorXd collectOnNodes(const std::vector<vec3d<double>>& nodesCloud,
                                   const operatorType inputOperatorType);

   private:
    double cc_;

    double getBasisValue(const vec3d<double>& nodeI, const vec3d<double>& nodeJ,
                         const operatorType inputOperatorType);

    // std::vector<double> NormVec_;
};

#endif
