
#ifndef GMMMQBASIS2D_H
#define GMMMQBASIS2D_H

#include<stdio.h>
#include<assert.h>
#include<iostream>
#include<math.h>
#include<vector>

using namespace std;

/*************************************************************************


*************************************************************************/
enum Operator_set
     {
        No_operation,
        Laplace,
        Partial_x,
        Partial_y
     };

class GMMMQBasis2D
{
    private:
        Operator_set Operator2D; 
        double CC;
        void Initialize();


    public:
        GMMMQBasis2D() : CC(1.0) {Initialize();};
        GMMMQBasis2D(double CCC) : CC(CCC) {Initialize();};
        ~GMMMQBasis2D() {}; 

        double SetOperator(const int OOperator2D,const vector<double>& VX1, const vector<double>& VX2);
};

void GMMMQBasis2D::Initialize()
{
    cout<<"MQ Basis 2D is chosen"<<endl;
    cout<<"Coefficient: "<<CC<<endl;       

}

double GMMMQBasis2D::SetOperator(const int OOperator2D, const vector<double>& VX1, const vector<double>& VX2)
{
    Operator2D = Operator_set(OOperator2D);

    double rx =  VX1[0] - VX2[0];
    double rz =  VX1[1] - VX2[1];
    double rs =  rx*rx+rz*rz;
    
    double temp;

    if (Operator2D == No_operation)
    {
        temp = sqrt(rs+CC*CC);
    }
    else if (Operator2D == Laplace)
    {
        double temp1;
        temp1 = sqrt(rs+CC*CC)*(rs+CC*CC);
        temp = (rs+2*CC*CC)/temp1;
    }
    else
    {
        cout << ">> Operator is not defined!" << endl;
        system("pause");
    }

    return temp;
};


#endif
