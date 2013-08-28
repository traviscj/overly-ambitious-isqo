
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
	
    {
        DenseAmplNlp dense_problem(problem_file);
        iSQOIterate penalty_iterate = dense_problem.initial();    
        iSQOQuadraticSubproblem dense_subproblem(dense_problem, penalty_iterate);
        SolveDenseQuadraticProgram dense_solver(dense_problem);
        iSQOStep dense_sol = dense_solver(dense_subproblem);
        std::cout << "dense sol: " << dense_sol << std::endl;   
    }
    {
        SparseAmplNlp sparse_problem(problem_file);
        iSQOIterate penalty_iterate = sparse_problem.initial();    
        iSQOSparseQuadraticSubproblem sparse_subproblem(sparse_problem, penalty_iterate);
        SolveSparseQuadraticProgram sparse_solver(sparse_problem);
        iSQOStep sparse_sol = sparse_solver(sparse_subproblem);
        std::cout << "sparse sol: " << sparse_sol << std::endl;   
    }
      
      // std::vector<double> example_solution(subproblem.num_qp_variables_);
      // example_.getPrimalSolution(&example_solution[0]);
      // std::vector<double> example_sparse_solution(sparse_subproblem.num_qp_variables_);
      // example_sparse.getPrimalSolution(&example_sparse_solution[0]);
      // 
      // // std::cout << "dense primal : " << example_solution << std::endl;
      // std::cout << "sparse primal: " << example_sparse_solution << std::endl;
}