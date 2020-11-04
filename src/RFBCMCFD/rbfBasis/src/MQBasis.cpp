#include "MQBasis.hpp"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

// for nuemann, any better idea?
// void MQBasis::setOperatorStatus(rbfOperatorType OOperatorStatus,
//                                   const std::vector<double>& NNormVec){
// inputOperatorType = OOperatorStatus;
// NormVec.resize( NNormVec.size() );
// copy( NNormVec, NormVec );
// };

MQBasis::MQBasis(std::shared_ptr<controlData> inControlData)
    : controlData_(inControlData)
{
    shapeParameter_ =
        controlData_->paramsDataAt({"solverControl", "shapeParameter"});

    dim_ = controlData_->paramsDataAt({"solverControl", "dimension"});
}

MQBasis::MQBasis(const double shapeParameter, const size_t dim)
    : shapeParameter_(shapeParameter), dim_(dim){};

double MQBasis::getBasisValue(const vec3d<double>& nodeI,
                              const vec3d<double>& nodeJ,
                              const rbfOperatorType inputOperatorType) const
{
    double rs = (nodeI - nodeJ).squaredNorm();
    vec3d<double> rr = nodeI - nodeJ;

    double temp;
    if (inputOperatorType == rbfOperatorType::CONSTANT)
    {
        temp = std::sqrt(rs + shapeParameter_ * shapeParameter_);
    }
    else if (inputOperatorType == rbfOperatorType::LAPLACE)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_) *
                (rs + shapeParameter_ * shapeParameter_);
        temp = ((dim_ - 1) * rs + dim_ * shapeParameter_ * shapeParameter_) /
               temp1;
    }
    else if (inputOperatorType == rbfOperatorType::DIVERGENCE)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = (rr[0] + rr[1] + rr[2]) / temp1;
    }
    else if (inputOperatorType == rbfOperatorType::PARTIAL_D1)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = rr[0] / temp1;
    }
    else if (inputOperatorType == rbfOperatorType::PARTIAL_D2)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = rr[1] / temp1;
    }
    else if (inputOperatorType == rbfOperatorType::PARTIAL_D3)
    {
        double temp1;
        temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
        temp = rr[2] / temp1;
    }
    // else if (inputOperatorType == rbfOperatorType::NEUMANN_OPERATOR)
    // {
    //     double temp1;
    //     temp1 = std::sqrt(rs + shapeParameter_ * shapeParameter_);
    //     temp = 0.0;

    //     for (size_t dim = 0; dim < dim; dim++)
    //     {
    //         temp += NormVec[dim] * rr[dim] / temp1;
    //     }
    // }
    else
    {
        std::cout << "Operator is not defined!" << std::endl;
        assert(false);
    }

    return temp;
}

Eigen::VectorXd MQBasis::collectOnNodes(
    const nodesCloud& cloud, const std::vector<vec3d<double>>& nodes,
    const rbfOperatorType operatorType) const
{
    std::cout << "cloooo" << std::endl;
    const size_t neighborNum = cloud.size_;
    Eigen::MatrixXd phi(neighborNum, neighborNum);

    for (size_t i = 0; i < neighborNum; i++)
    {
        for (size_t j = 0; j < neighborNum; j++)
        {
            phi(j, i) /* transposed */ =
                getBasisValue(nodes[cloud.ids_[i]], nodes[cloud.ids_[j]],
                              rbfOperatorType::CONSTANT);
        }
    }

    Eigen::VectorXd phiL(neighborNum);
    for (size_t j = 0; j < neighborNum; j++)
    {
        phiL(j) = getBasisValue(nodes[cloud.id0_], nodes[cloud.ids_[j]],
                                operatorType);
    }

    Eigen::VectorXd x(neighborNum);
    x = phi.ldlt().solve(phiL);
    std::cout << "cloooo===============" << std::endl;

    return x;
}
