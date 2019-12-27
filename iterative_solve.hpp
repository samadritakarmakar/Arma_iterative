#ifndef ITERATIVE_SOLVE_HPP
#define ITERATIVE_SOLVE_HPP
#include "armadillo"
#include "lis.h"
#include "lis_matrix.h"
#include "lis_vector.h"
#include <string>
using namespace arma;

enum class precond {none, jacobi, ilu, ssor, hybrid, is, sainv, saamg, iluc, ilut};

enum class solver_type {CG, BiCG, CGS, BiCGSTAB, BiCGSTAB_l, GPBiCG, TFQMR, Orthmin, GMRES, Jacobi, GaussSeidel, SOR,
                        BiCG_Safe, CR, BiCR, CRS, BiCRSTAB, GPBiCR, BiCR_Safe, FGMRES, IDRs, IDRl, MINRES, COCG, COCR};

std::string GetPreconditioner(precond preconditioner=precond::jacobi)
{
    std::string optionb4precon;
    optionb4precon=" -p ";
    switch (preconditioner)
    {
    case precond::none:
        optionb4precon+="none";
        break;
    case precond::jacobi:
        optionb4precon+="jacobi";
        break;
    case precond::ilu:
        optionb4precon+="ilu";
        break;
    case precond::ssor:
        optionb4precon+="ssor";
        break;
     case precond::hybrid:
        optionb4precon+="hybrid";
        break;
     case precond::is:
        optionb4precon+="is";
        break;
    case precond::sainv:
       optionb4precon+="sainv";
       break;
    case precond::saamg:
       optionb4precon+="saamg";
       break;
     case precond::iluc:
        optionb4precon+="iluc";
        break;
    case precond::ilut:
        optionb4precon+="ilut";
        break;
    default:
        cout<<"Precoditioner not found!\n";
        throw;
        break;
    }
    return (optionb4precon);
}

std::string GetSolver(solver_type solver=solver_type::BiCG)
{
    std::string option4solver="-i ";
    switch (solver)
    {
    case solver_type::CG:
        option4solver+="cg";
        break;
    case solver_type::BiCG:
        option4solver+="bicg";
        break;
    case solver_type::CGS:
        option4solver+="cgs";
        break;
    case solver_type::BiCGSTAB:
        option4solver+="bicgstab";
        break;
    case solver_type::BiCGSTAB_l:
        option4solver+="bicgstabl";
        break;
    case solver_type::GPBiCG:
        option4solver+="gpbicg";
        break;
    case solver_type::TFQMR:
        option4solver+="tfqmr";
        break;
    case solver_type::Orthmin:
        option4solver+="orthomin";
        break;
    case solver_type::GMRES:
        option4solver+="gmres";
        break;
    case solver_type::Jacobi:
        option4solver+="jacobi";
        break;
    case solver_type::GaussSeidel:
        option4solver+="gs";
        break;
    case solver_type::SOR:
        option4solver+="sor";
        break;
    case solver_type::BiCG_Safe:
        option4solver+="bicgsafe";
        break;
    case solver_type::CR:
        option4solver+="cr";
        break;
    case solver_type::BiCR:
        option4solver+="bicr";
        break;
    case solver_type::CRS:
        option4solver+="crs";
        break;
    case solver_type::BiCR_Safe:
        option4solver+="bicrsafe";
        break;
    case solver_type::FGMRES:
        option4solver+="fgmres";
        break;
    case solver_type::IDRs:
        option4solver+="idrs";
        break;
    case solver_type::IDRl:
        option4solver+="idrl";
        break;
    case solver_type::MINRES:
        option4solver+="minres";
        break;
    case solver_type::COCG:
        option4solver+="cocg";
        break;
    case solver_type::COCR:
        option4solver+="cocr";
        break;
    default:
        break;
    }
    return option4solver;
}
template<class vector>
void iterative_solve(mat& x, sp_mat& A, vector&b, std::string solver_tolerance="1e-12",
                     solver_type solver_typ=solver_type::GMRES, precond preconditioner=precond::ilu)
{
    A.sync();
    int argc=0; char** argv=NULL;
    lis_initialize(&argc, &argv);
    LIS_SCALAR *value;
    LIS_MATRIX A1;
    LIS_SOLVER solver;
    LIS_INT nnz;
    //LIS_INT ptr[A.n_cols+1],index[A.n_nonzero];
    LIS_INT *ptr,*index;
    nnz=A.n_nonzero;
    value= const_cast<LIS_SCALAR *>(reinterpret_cast <const LIS_SCALAR *>(A.values));
    ptr = const_cast<LIS_INT *>(reinterpret_cast<const LIS_INT *>(A.col_ptrs));
    index = const_cast<LIS_INT *>(reinterpret_cast<const LIS_INT *>(A.row_indices));
    lis_matrix_create(0,&A1);
    lis_matrix_set_size(A1,0,(int)A.n_rows);
    lis_matrix_set_csc(nnz,ptr,index,value,A1);
    lis_matrix_assemble(A1);
    LIS_SCALAR *value_b, *value_x;
    LIS_VECTOR b1,x1;
    lis_vector_duplicate(A1,&b1);
    lis_vector_duplicate(A1,&x1);
    if(x.n_rows <A.n_rows || x.n_cols!=1)
    {
        x.set_size(A.n_rows, 1);
        x.zeros();
    }

    value_b= const_cast<LIS_SCALAR *>(reinterpret_cast <const LIS_SCALAR *>(b.memptr()));
    value_x= const_cast<LIS_SCALAR *>(reinterpret_cast <const LIS_SCALAR *>(x.memptr()));
    lis_vector_set(b1, value_b);
    lis_vector_set(x1, value_x);
    lis_solver_create(&solver);
    std::string solver_option1_string=GetSolver(solver_typ)+GetPreconditioner(preconditioner);
    char *solver_option1=const_cast<char *>(solver_option1_string.c_str());
    lis_solver_set_option(solver_option1,solver);
    std::string tolerance="-tol "+solver_tolerance+" -maxiter "+std::to_string(A.n_rows) +" -print out";
    //cout<<tolerance<<"\n";
    char *solver_option2=const_cast<char *>(tolerance.c_str());
    lis_solver_set_option(solver_option2,solver);
    int success =lis_solve(A1,b1,x1,solver);
    if(success) //if unsuccessful, it returns a number other thatn zero, and hence the program stops.
    {
        std::cout<<"success= "<<success<<"\n";
        std::cout<<"Could not solve! Please try a iterative different solver config or use a direct solver.\n";
        throw;
    }
    //bool copy_aux_mem = false;
    //bool strict = false;
    //x=mat(x1->value, A.n_rows, 1, copy_aux_mem, strict);
    lis_finalize();
}

template<class vector>
vector iterative_solve(sp_mat& A, vector&b, std::string solver_tolerance="1e-12",
                    solver_type solver_typ=solver_type::GMRES, precond preconditioner=precond::ilu)
{
    vector x;
    iterative_solve(x, A, b, solver_tolerance, solver_typ, preconditioner);
    return x;
}

/*void iterative_solve(mat& x, sp_mat& A, mat&b, std::string solver_tolerance="1e-12",
                    solver_type solver_typ=solver_type::GMRES, precond preconditioner=precond::ilu)
{
    bool copy_aux_mem = false;
    bool strict1 = true;
    bool strict2 = false;
    vec x1;
    vec b1;
    b1=mat(b.memptr(), A.n_rows, 1, copy_aux_mem, strict1);
    iterative_solve(x1, A, b1, solver_tolerance, solver_typ, preconditioner);
    //cout<<x1;
    x=mat(x1.memptr(), A.n_rows, 1, copy_aux_mem, strict1);
}

mat iterative_solve(sp_mat& A, mat&b, std::string solver_tolerance="1e-12",
                    solver_type solver_typ=solver_type::GMRES, precond preconditioner=precond::ilu)
{
    bool copy_aux_mem = false;
    bool strict = false;
    mat x;
    iterative_solve(x, A, b, solver_tolerance, solver_typ, preconditioner);
    //mat x=mat(x1.memptr(), A.n_rows, 1, copy_aux_mem, strict);
    return x;
}*/

#endif // ITERATIVE_SOLVE_HPP
