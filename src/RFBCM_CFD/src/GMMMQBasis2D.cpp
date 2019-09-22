#include "GMMMQBasis2D.h"
#include "gmm.h"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

using namespace std;
using namespace gmm;

void GMMMQBasis2D::SetOperatorStatus( OperatorType OOperatorStatus ) {
    OperatorStatus = OOperatorStatus;
};

// for nuemann, any better idea?
void GMMMQBasis2D::SetOperatorStatus( OperatorType            OOperatorStatus,
                                      const vector< double >& NNormVec ) {
    OperatorStatus = OOperatorStatus;
    NormVec.resize( NNormVec.size() );
    copy( NNormVec, NormVec );
};

double GMMMQBasis2D::GetBasisValue( const vector< double >& VX1, const vector< double >& VX2,
                                    int flag ) {
    int DIMENSION = VX1.size();

    vector< double > rr( DIMENSION );

    add( VX1, scaled( VX2, -1.0 ), rr );

    double rs = vect_sp( rr, rr );

    double temp;

    //
    if ( OperatorStatus == OperatorType::IdentityOperation || flag != 0 ) {
        temp = sqrt( rs + CC * CC );
    }
    else if ( OperatorStatus == OperatorType::Laplace ) {
        double temp1;
        temp1 = sqrt( rs + CC * CC ) * ( rs + CC * CC );
        temp  = ( rs + 2 * CC * CC ) / temp1;
    }
    else if ( OperatorStatus == OperatorType::Partial_D1 ) {
        double temp1;
        temp1 = sqrt( rs + CC * CC );
        temp  = rr[ 0 ] / temp1;
    }
    else if ( OperatorStatus == OperatorType::Partial_D2 ) {
        double temp1;
        temp1 = sqrt( rs + CC * CC );
        temp  = rr[ 1 ] / temp1;
    }
    else if ( OperatorStatus == OperatorType::NEUMANN_OPERATOR ) {
        double temp1;
        temp1 = sqrt( rs + CC * CC );
        temp  = 0.0;

        for ( int dim = 0; dim < DIMENSION; dim++ ) {
            temp += NormVec[ dim ] * rr[ dim ] / temp1;
        }
    }
    else {
        cout << "Operator is not defined!" << endl;
        assert( false );
    }

    return temp;
};
