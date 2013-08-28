
// standard libraries:
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

// iSQO headers:
#include "isqo.hh"

int main(int argc, char **argv) {
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs067.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
	AmplNlp problem(problem_file);
	iSQOIterate penalty_iterate = problem.initial();    
    for (int ind=0; ind < penalty_iterate.num_dual_eq_; ++ind) penalty_iterate.dual_eq_values_[ind] = 0.0;
    for (int ind=0; ind < penalty_iterate.num_dual_ieq_; ++ind) penalty_iterate.dual_ieq_values_[ind] = 0.0;
	
    std::cout << "objective: " << problem.objective(penalty_iterate) << std::endl;
	std::cout << "eq constraints: "<< problem.constraints_equality(penalty_iterate) << std::endl;
	
    sparse_matrix sparse_eq_jacobian = problem.constraints_equality_jacobian_sparse(penalty_iterate);
    std::cout << "eq jacobian: " << sparse_eq_jacobian << std::endl;
    matrix eq_jacobian = problem.constraints_equality_jacobian(penalty_iterate);
    std::cout << "eq jacobian: " << eq_jacobian << std::endl;
    
    
    std::cout << "ieq constraints: "<< problem.constraints_inequality(penalty_iterate) << std::endl;
    sparse_matrix sparse_ieq_jacobian = problem.constraints_inequality_jacobian_sparse(penalty_iterate);
    std::cout << "ieq jacobian: " << sparse_ieq_jacobian << std::endl;
    matrix ieq_jacobian = problem.constraints_inequality_jacobian(penalty_iterate);
    std::cout << "ieq jacobian: " << ieq_jacobian << std::endl;
    
    sparse_matrix hessian = problem.lagrangian_hessian_sparse(penalty_iterate);
    std::cout << "hessian: " << hessian << std::endl;
    
    
    
    // double *hess_full = qpOASES_hessian.full();
    // for (size_t ind=0; ind < hessian.num_rows() * hessian.num_columns(); ++ind) {
    //     std::cout << "H[" << ind << "]: " << hess_full[ind] << std::endl;
    // }
    
    iSQOQuadraticSubproblem subproblem(problem, penalty_iterate);
    std::cout << "================" << std::endl;
    std::cout << "hessian: " << subproblem.hessian_sparse_ << std::endl;
    std::cout << "================" << std::endl;
    
    qpOASES::SQProblem example_(subproblem.num_qp_variables_, subproblem.num_qp_constraints_, qpOASES::HST_SEMIDEF);
    int nWSR=10000;
	int ret = example_.init( &subproblem.hessian_.data_[0],
						 &subproblem.gradient_[0],
						 &subproblem.jacobian_.data_[0],
						 &subproblem.lower_bound_[0],
						 &subproblem.upper_bound_[0],
						 &subproblem.jacobian_lower_bound_[0],
						 &subproblem.jacobian_upper_bound_[0],
						 nWSR,
						 0);
     std::cout << "ret: " << ret << ", " << nWSR << std::endl;
    
     std::cout << "sparse jacobian: " << subproblem.jacobian_sparse_ << std::endl;
     qpOASES::SparseMatrix qpOASES_jacobian(subproblem.num_qp_constraints_, subproblem.num_qp_variables_, &subproblem.jacobian_sparse_.row_indices_[0], &subproblem.jacobian_sparse_.col_starts_[0], &subproblem.jacobian_sparse_.vals_[0]);
 	double* qpOASES_jacobian_full = qpOASES_jacobian.full();
     for (size_t ind=0; ind < subproblem.num_qp_constraints_*subproblem.num_qp_variables_; ++ind) {
         // std::cout << "qpOASES_jacobian_full[" << ind << "]: " << [ind] << std::endl;
         // if (subproblem.jacobian_.data_[ind] != qpOASES_jacobian_full[ind])
          {
             std::cout  << "qpOASES_jacobian_full[" << ind 
                        << "]: dense: " << subproblem.jacobian_.data_[ind] 
                        << "; sparse: " << qpOASES_jacobian_full[ind] << std::endl;
         }
     }
     qpOASES::SymSparseMat qpOASES_hessian(subproblem.num_qp_variables_, subproblem.num_qp_variables_, &subproblem.hessian_sparse_.row_indices_[0], &subproblem.hessian_sparse_.col_starts_[0], &subproblem.hessian_sparse_.vals_[0]);
     // qpOASES_hessian.createDiagInfo();
 	double* qpOASES_hessian_full = qpOASES_hessian.full();
     for (size_t ind=0; ind < subproblem.num_qp_variables_*subproblem.num_qp_variables_; ++ind) {
         // if (subproblem.hessian_.data_[ind] != qpOASES_hessian_full[ind]) 
         {
             std::cout  << "qpOASES_hessian_full[" << ind 
                        << "]: dense: " << subproblem.hessian_.data_[ind] 
                        << "; sparse: " << qpOASES_hessian_full[ind] << std::endl;
         }
     }
     qpOASES::SQProblem example_sparse(subproblem.num_qp_variables_, subproblem.num_qp_constraints_, qpOASES::HST_SEMIDEF);
     int nWSR2=10000;
 	int ret2 = example_sparse.init( &qpOASES_hessian,
 						 &subproblem.gradient_[0],
 						 &qpOASES_jacobian,
 						 &subproblem.lower_bound_[0],
 						 &subproblem.upper_bound_[0],
 						 &subproblem.jacobian_lower_bound_[0],
 						 &subproblem.jacobian_upper_bound_[0],
 						 nWSR2,
 						 0);
      std::cout << "ret2: " << ret2 << ", " << nWSR2 << std::endl;
      
      std::vector<double> example_solution(subproblem.num_qp_variables_);
      example_.getPrimalSolution(&example_solution[0]);
      std::vector<double> example_sparse_solution(subproblem.num_qp_variables_);
      example_sparse.getPrimalSolution(&example_sparse_solution[0]);
      
      std::cout << "dense primal : " << example_solution << std::endl;
      std::cout << "sparse primal: " << example_sparse_solution << std::endl;
}