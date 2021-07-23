#ifndef MQBasis_HPP
#define MQBasis_HPP

#include "MeshData.hpp"
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
    MQBasis(ControlData* controlData, MeshData* meshData);
    MQBasis(const double shapeParameter, const size_t dim);
    ~MQBasis(){};

    Eigen::VectorXd collectOnNodes(const size_t& nodeID,
                                   const rbfOperatorType& operatorType) const;

   private:
    ControlData* controlData_;
    MeshData* meshData_;
    double shapeParameter_;
    size_t dim_;

    double getBasisValue(const size_t& i, const size_t& j,
                         const rbfOperatorType& operatorType) const;

    // std::vector<double> NormVec_;
};

#endif
