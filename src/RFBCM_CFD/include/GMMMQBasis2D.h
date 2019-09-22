#ifndef GMMMQBASIS2D_H
#define GMMMQBASIS2D_H

#include <vector>

using namespace std;

/*************************************************************************


*************************************************************************/
class GMMMQBasis2D {

public:
    enum class OperatorType {
        IdentityOperation,
        Laplace,
        Partial_D1,
        Partial_D2,
        NEUMANN_OPERATOR
    };

    GMMMQBasis2D( double CCC = 1.0 );
    ~GMMMQBasis2D(){};

    void SetOperatorStatus( OperatorType OOperatorStatus );
    void SetOperatorStatus( OperatorType OOperatorStatus, const vector< double >& NNormVec );

    double GetBasisValue( const vector< double >& VX1, const vector< double >& VX2, int flag = 0 );

private:
    double           CC;
    OperatorType     OperatorStatus;
    vector< double > NormVec;
};

#endif
