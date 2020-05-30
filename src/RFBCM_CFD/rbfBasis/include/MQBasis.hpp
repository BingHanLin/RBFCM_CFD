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
    MQBasis(std::shared_ptr<controlData> inControlData);
    MQBasis(const double shapeParameter, const int dim);
    ~MQBasis(){};

    Eigen::VectorXd collectOnNodes(
        const std::vector<vec3d<double>>& nodesCloud,
        const rbfOperatorType inputOperatorType) const;

   private:
    std::shared_ptr<controlData> controlData_;
    double shapeParameter_;
    int dim_;

    double getBasisValue(const vec3d<double>& nodeI, const vec3d<double>& nodeJ,
                         const rbfOperatorType inputOperatorType) const;

    // std::vector<double> NormVec_;
};

#endif
