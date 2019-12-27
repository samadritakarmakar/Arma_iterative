# Arma Iterative is an attempt to extend the capabilites of [Armadillo](http://arma.sourceforge.net/) library by linking it with [lis](https://www.ssisc.org/lis/) Solver Library with minute modifications.  
The lis solver library provides a wide range of iterative solvers and supports compressed sparse column (CSC) format. This makses it easy to integrate with the armadillo library. The only contribution from me in this library is the iterative_solve.hpp file. Rest is left to lis to take care of.  
You can use this library by calling iterative_solve function in the following way.  

       iterative_solve(x, A, b, "1e-12", solver_type::BiCGSTAB, precond::hybrid);
       
The vector x is the solution vector of the equation Ax=b, where A is the sparse matrix, b is the right hand side vector. "1e-12" is the tolerance used. It has to be passed in form of a string.  
The solver options available are CG, BiCG, CGS, BiCGSTAB, BiCGSTAB_l, GPBiCG, TFQMR, Orthmin, GMRES, Jacobi, GaussSeidel, SOR,
                        BiCG_Safe, CR, BiCR, CRS, BiCRSTAB, GPBiCR, BiCR_Safe, FGMRES, IDRs, IDRl, MINRES, COCG, COCR.  
The precondtioners availalble are  jacobi, ilu, ssor, hybrid, is, sainv, saamg, iluc, ilut.  
For more details, please see [lis](lis/doc/lis-ug-en.pdf) documentation.  
I am making this small source code availaible under Apache license. Armadillo and lis are licensed under Apache and BSD licenses.  
