#include "MQBasis.hpp"
#include "messages.hpp"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

// for nuemann, any better idea?
// void MQBasis::setOperatorStatus(rbfOperatorType OOperatorStatus,
//                                   const std::vector<double>& NNormVec){
// operatorType = OOperatorStatus;
// NormVec.resize( NNormVec.size() );
// copy( NNormVec, NormVec );
// };

MQBasis::MQBasis(ControlData* controlData, MeshData* meshData)
    : controlData_(controlData), meshData_(meshData)
{
    shapeParameter_ =
        controlData_->paramsDataAt({"solverControl", "shapeParameter"});

    dim_ = controlData_->paramsDataAt({"solverControl", "dimension"});
}

MQBasis::MQBasis(const double shapeParameter, const size_t dim)
    : shapeParameter_(shapeParameter), dim_(dim){};

double MQBasis::getBasisValue(const size_t& i, const size_t& j,
                              const rbfOperatorType& operatorType) const
{
    const auto nodeI = meshData_->nodeByID(i);
    const auto nodeJ = meshData_->nodeByID(j);

    double rs = (nodeI - nodeJ).squaredNorm();
    vec3d<double> rr = nodeI - nodeJ;

    double temp;
    if (operatorType == rbfOperatorType::CONSTANT)
    {
        temp = std::sqrt(rs + shapeParameter_ * shapeParameter_);
    }
    else if (operatorType == rbfOperatorType::LAPLACE)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_) *
                (rs + shapeParameter_ * shapeParameter_);
        temp = ((dim_ - 1) * rs + dim_ * shapeParameter_ * shapeParameter_) /
               temp1;
    }
    else if (operatorType == rbfOperatorType::DIVERGENCE)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = (rr[0] + rr[1] + rr[2]) / temp1;
    }
    else if (operatorType == rbfOperatorType::PARTIAL_D1)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = rr[0] / temp1;
    }
    else if (operatorType == rbfOperatorType::PARTIAL_D2)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = rr[1] / temp1;
    }
    else if (operatorType == rbfOperatorType::PARTIAL_D3)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = rr[2] / temp1;
    }
    else if (operatorType == rbfOperatorType::NEUMANN)
    {
        const auto norm = meshData_->normalByID(i);

        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = 0.0;

        for (size_t d = 0; d < dim_; d++)
        {
            temp += norm[d] * rr[d] / temp1;
        }
    }
    else
    {
        ASSERT("Operator is not defined!");
    }

    return temp;
}

Eigen::VectorXd MQBasis::collectOnNodes(
    const size_t& nodeID, const rbfOperatorType& operatorType) const
{
    const auto cloud = meshData_->cloudByID(nodeID);

    const size_t neighborNum = cloud.size_;
    Eigen::MatrixXd phi(neighborNum, neighborNum);

    for (size_t i = 0; i < neighborNum; i++)
    {
        for (size_t j = 0; j < neighborNum; j++)
        {
            phi(j, i) /* transposed */ = getBasisValue(
                cloud.ids_[i], cloud.ids_[j], rbfOperatorType::CONSTANT);
        }
    }

    Eigen::VectorXd phiL(neighborNum);
    for (size_t j = 0; j < neighborNum; j++)
    {
        phiL(j) = getBasisValue(cloud.id0_, cloud.ids_[j], operatorType);
    }

    Eigen::VectorXd x(neighborNum);
    x = phi.ldlt().solve(phiL);

    return x;
}