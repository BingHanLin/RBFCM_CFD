#ifndef MQBASIS2D_HPP
#define MQBASIS2D_HPP

#include <Eigen/Dense>
#include <vector>

/*************************************************************************


*************************************************************************/
class MQBasis2D
{
   public:
    enum operatorType
    {
        CONSTANT,
        LAPLACE,
        PARTIAL_D1,
        PARTIAL_D2,
        NEUMANN_OPERATOR
    };

    MQBasis2D(double CCC = 1.0) : cc_(CCC){};

    ~MQBasis2D(){};

    Eigen::VectorXd collectOnNodes(
        const std::vector<std::vector<double>>& nodesCloud,
        const operatorType inputOperatorType);

   private:
    double cc_;

    double getBasisValue(const std::vector<double>& nodeI,
                         const std::vector<double>& nodeJ,
                         const operatorType inputOperatorType);

    // std::vector<double> NormVec_;
};

#endif
