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
    MQBasis(const double shapeParameter, const size_t dim);
    ~MQBasis(){};

    Eigen::VectorXd collectOnNodes(const nodesCloud& cloud,
                                   const std::vector<vec3d<double>>& nodes,
                                   const rbfOperatorType operatorType) const;

    Eigen::VectorXd collectOnNodes(const nodesCloud& cloud,
                                   const std::vector<vec3d<double>>& nodes,
                                   const vec3d<double>& norm,
                                   const rbfOperatorType operatorType) const;

   private:
    std::shared_ptr<controlData> controlData_;
    double shapeParameter_;
    size_t dim_;

    double getBasisValue(const vec3d<double>& nodeI, const vec3d<double>& nodeJ,
                         const vec3d<double>& norm,
                         const rbfOperatorType operatorType) const;

    // std::vector<double> NormVec_;
};

#endif
