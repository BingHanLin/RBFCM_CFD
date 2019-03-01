
#ifndef GMMRECTANGLE_H
#define GMMRECTANGLE_H

#include<stdio.h>
#include<assert.h>
#include<iostream>
#include<math.h>
#include<vector>

using namespace std;

/*************************************************************************


*************************************************************************/

class GMMRECTANGLE
{
    private:
        int M, N;
        double LX, LY;
        double CX, CY;

    public:
        int Nint, Nbou, Nall, Dim;
        GMMRECTANGLE (int mm, int nn, double LLx, double LLy);
        ~GMMRECTANGLE() {};

        vector<double> GetInteriorNode(int i) const;
        vector<double> GetBoundaryNode(int i) const;
        vector<double> GetAllNode(int i) const;        
        vector<double> GetOutNormal(int i) const;

};

GMMRECTANGLE::GMMRECTANGLE(int mm, int nn, double LLx, double LLy)
{
    M = mm;
    N = nn;
    LX = LLx;
    LY = LLy;
    CX = 0;
    CY = 0;
    Dim = 2;  // dimension 

    Nint = M*N;
    Nbou = 2*M+2*N;
    Nall = Nint + Nbou;

    cout<<"Rectangular mesh is created"<<endl;
    cout<<"Nint: "<<Nint<<", Nbou: "<<Nbou<<", Nall: "<<Nall<<endl;
}

vector<double> GMMRECTANGLE::GetInteriorNode(int i) const
{
    assert(i>=1 && i<=M*N);
    vector<double> temp(2);

    int I = (i-1)/M;
    int J = (i-1)%M;

    temp[0] = ((double)(J)+0.5)*LX/(double)(M);
    temp[1] = ((double)(I)+0.5)*LY/(double)(N);

    temp[0] = temp[0]-LX/2.+CX;
    temp[1] = temp[1]-LY/2.+CY;

    return temp;
}

vector<double> GMMRECTANGLE::GetBoundaryNode(int i) const
{
    assert(i>=1 && i<=2*(N+M));
    vector<double> temp(2);

    if(i<=M)
    {
        temp[0]=(double)(i-1)*LX/(double)(M)+.5*LX/(double)(M);
        temp[1]=LY;
    }
    else if(i<=M+N)
    {
        temp[0]=LX;
        temp[1]=LY-.5*LY/(double)(N)-(double)(i-M-1)*LY/(double)(N);
    }
    else if(i<=M+N+M)
    {
        temp[0]=LX-.5*LX/(double)(M)-(double)(i-M-N-1)*LX/(double)(M);
        temp[1]=0.;
    }
    else
    {
        temp[0]=0.;
        temp[1]=(double)(i-M-N-M-1)*LY/(double)(N)+.5*LY/(double)(N);
    }

    temp[0]=temp[0]-LX/2.+CX;
    temp[1]=temp[1]-LY/2.+CY;

    return temp;
}

vector<double> GMMRECTANGLE::GetAllNode(int i) const
{
    assert(i>=1 && i<=M*N+2*(N+M));
    vector<double> temp(2);

    if(i <= M*N)
        return GetInteriorNode(i);
    else
        return GetBoundaryNode(i-M*N);
}


vector<double> GMMRECTANGLE::GetOutNormal(int i) const
{
    assert( i>=1 && i<=2*(N+M) );
    vector<double> temp(2);

    if(i<=M)
    {
        temp[0] = 0.0;
        temp[1] = 1.0;
    }
    else if(i<=M+N)
    {
        temp[0] = 1.0;
        temp[1] = 0.0;
    }
    else if(i<=M+N+M)
    {
        temp[0]= 0.0;
        temp[1]= -1.0;
    }
    else
    {
        temp[0]= -1.0;
        temp[1]= 0.0;
    }

    return temp;
}



#endif
