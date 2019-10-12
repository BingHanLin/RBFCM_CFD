#ifndef MQBASIS2D_HPP
#define MQBASIS2D_HPP

#include <vector>

using namespace std;

/*************************************************************************


*************************************************************************/
class MQBasis2D {

public:
    enum operatorType { IdentityOperation, Laplace, Partial_D1, Partial_D2, NEUMANN_OPERATOR };

    MQBasis2D( double CCC = 1.0 ) : cc_( CCC ){};
    ~MQBasis2D(){};

    void setOperatorStatus( operatorType OOperatorStatus );
    void setOperatorStatus( operatorType OOperatorStatus, const std::vector< double >& NNormVec );

    double getBasisValue( const std::vector< double >& VX1, const std::vector< double >& VX2,
                          int flag = 0 );
    double cc_;

private:
    operatorType          operatorStatus_;
    std::vector< double > NormVec_;
};

#endif
