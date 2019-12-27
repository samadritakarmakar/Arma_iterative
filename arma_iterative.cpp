#ifndef ARMA_ITERATIVE_HPP
#define ARMA_ITERATIVE_HPP
#include "iterative_solve.hpp"

using namespace arma;
int main()
{
    sp_mat A(5,5);
    mat b(5,1);
    for(double i=1.;i<=A.n_rows; i++)
    {
        b(i-1,0)=i;
    }
    mat x;
    sp_mat M(5,5);
    /*A= 4     1      2    0.5   2.0
         1     0.5    0     0     0
         2     0      3.0   0     2
         0.5   0      0   .625    0
         2.0   0      2     0     16*/

    A(0,1)=1;
    A(0,2)=2;
    A(0,3)=2;
    A(0,3)=0.5;
    A(0,4)=2.0;
    A+=A.t();

    A(0,0)=4;
    A(1,1)=0.5;
    A(1,2)=0.;
    A(2,1)=0.;
    A(2,2)=3.0;
    A(3,3)=0.625;
    A(4,4)=16;
    A(4,2)=2;
    A(2,4)=2;

    iterative_solve(x, A, b);
    vec x1;
    spsolve(x1,A,b);
    cout<<x1<<"\n";
    cout<<x<<"\n";
    return 0;
}

#endif // ARMA_ITERATIVE_HPP
