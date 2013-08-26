
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
	
    std::cout << "objective: " << problem.objective(penalty_iterate) << std::endl;
	std::cout << "constraints: "<< problem.constraints_equality(penalty_iterate) << std::endl;
	
    qpOASES::SparseMatrix sparse_jacobian = problem.constraints_equality_jacobian_sparse(penalty_iterate);
    // std::cout << "eq jacobian: " << ) << std::endl;
    
}