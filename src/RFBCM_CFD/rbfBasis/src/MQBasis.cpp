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
        controlData_->paramsDataAt({"solverConstrol", "shapeParameter"});
}

MQBasis::MQBasis(const double shapeParameter)
    : shapeParameter_(shapeParameter){};

double MQBasis::getBasisValue(const vec3d<double>& nodeI,
                              const vec3d<double>& nodeJ,
                              const rbfOperatorType inputOperatorType)
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
        temp = (rs + 2 * shapeParameter_ * shapeParameter_) / temp1;
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

    //     for (int dim = 0; dim < dim; dim++)
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
    const std::vector<vec3d<double>>& nodesCloud,
    const rbfOperatorType inputOperatorType)
{
    const int neighborNum = nodesCloud.size();
    Eigen::MatrixXd phi(neighborNum, neighborNum);

    for (int i = 0; i < neighborNum; i++)
    {
        for (int j = 0; j < neighborNum; j++)
        {
            phi(j, i) /* transposed */ = getBasisValue(
                nodesCloud[i], nodesCloud[j], rbfOperatorType::CONSTANT);
        }
    }

    Eigen::VectorXd phiL(neighborNum);
    for (int j = 0; j < neighborNum; j++)
    {
        phiL(j) =
            getBasisValue(nodesCloud[0], nodesCloud[j], inputOperatorType);
    }

    Eigen::VectorXd x(neighborNum);
    x = phi.ldlt().solve(phiL);

    return x;
}
