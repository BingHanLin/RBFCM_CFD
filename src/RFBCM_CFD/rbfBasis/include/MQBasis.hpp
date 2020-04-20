#ifndef MQBasis_HPP
#define MQBasis_HPP

#include "controlData.hpp"
#include "dataStructure.hpp"
#include "enumMap.hpp"

#include <Eigen/Dense>
#include <vector>

/*************************************************************************


*************************************************************************/
class MQBasis
{
   public:
    MQBasis();

    ~MQBasis(){};

    Eigen::VectorXd collectOnNodes(const std::vector<vec3d<double>>& nodesCloud,
                                   const rbfOperatorType inputOperatorType);

   private:
    controlData* controlData_;
    double shapeParameter_;

    double getBasisValue(const vec3d<double>& nodeI, const vec3d<double>& nodeJ,
                         const rbfOperatorType inputOperatorType);

    // std::vector<double> NormVec_;
};

#endif
