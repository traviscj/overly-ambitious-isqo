// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

// standard libraries:
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

// iSQO headers:
#include "isqo.hh"

int main() {
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs009.nl");
	AmplNlp problem(problem_file);
	iSQOIterate penalty_iterate = problem.initial();
	matrix je = problem.constraints_equality_jacobian(penalty_iterate);
	std::vector<double> x(penalty_iterate.primal_values_);
	std::vector<double> y(penalty_iterate.dual_eq_values_);
	std::cout << "jac: "<< je << std::endl;
	std::cout << "x: " << x << std::endl;
	std::cout << je.multiply(x) << std::endl;
	
	std::cout << "y: " << y << std::endl;
	std::cout << je.multiply_transpose(y) << std::endl;
	
}