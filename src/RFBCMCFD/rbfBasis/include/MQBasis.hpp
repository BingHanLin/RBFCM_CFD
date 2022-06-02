#ifndef MQBasis_HPP
#define MQBasis_HPP

#include "controlData.hpp"
#include "meshData.hpp"

#include "dataStructure.hpp"
#include "enumMap.hpp"

#include <Eigen/Dense>
#include <vector>

/*************************************************************************


*************************************************************************/
class MQBasis
{
   public:
    MQBasis(const double& shapeParameter, const size_t& didimension,
            MeshData* meshData);
    ~MQBasis(){};

    Eigen::VectorXd collectOnNodes(const size_t& nodeID,
                                   const rbfOperatorType& operatorType) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    double shapeParameter_;
    size_t dimension_;

    double getBasisValue(const size_t& i, const size_t& j,
                         const rbfOperatorType& operatorType) const;

    // std::vector<double> NormVec_;
};

#endif
