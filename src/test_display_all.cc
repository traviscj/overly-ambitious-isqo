
// standard libraries:
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

// iSQO headers:
#include "isqo.hh"

int main() {
	AmplNlp ampl_nlp_object("/Users/traviscj/optimization/cute_nl_nopresolve/hs014.nl");
	
	
	Hs014 problem;
	iSQOIterate iterate(2,1,1);
	iterate.penalty_parameter_ = 1.0e-1;
	iterate.primal_values_[0] = 2.0;
	iterate.primal_values_[1] = 2.0;
	iterate.dual_eq_values_[0] = -1.0;
	iterate.dual_ieq_values_[0] = 1000;
	// iterate.dual_ieq_values_[0] = ;
	// iterate.dual_ieq_values_[0] = 0;
	// iterate.penalty_parameter_ = 1.0e-1;
	// iterate.primal_values_[0] = .822875656;
	// iterate.primal_values_[1] = .911437828;
	// iterate.dual_eq_values_[0] = 1.846589027861980e+00;
	// iterate.dual_ieq_values_[0] = 1.594493103554523e+00;

	std::cout << "Zeroth Order:" << std::endl;
	std::cout << "===[Objective]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.objective(iterate) << std::endl;
	std::cout << "Custom : " << problem.objective(iterate) << std::endl;
	
	std::cout << "===[Equality Constraints]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_equality(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_equality(iterate) << std::endl;
	
	std::cout << "===[Inequality Constraints]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_inequality(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_inequality(iterate) << std::endl;
	
	std::cout << "First Order:" << std::endl;
	std::cout << "===[Objective Gradient]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.objective_gradient(iterate) << std::endl;
	std::cout << "Custom : " << problem.objective_gradient(iterate) << std::endl;
	
	std::cout << "===[Equality Jacobian]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_equality_jacobian(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_equality_jacobian(iterate) << std::endl;
	
	
	std::cout << "===[Inequality Jacobian]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.constraints_inequality_jacobian(iterate) << std::endl;
	std::cout << "Custom : " << problem.constraints_inequality_jacobian(iterate) << std::endl;
	
	std::cout << "Second Order:" << std::endl;
	std::cout << "===[Lagrangian Hessian]==============" << std::endl;
	std::cout << "  AMPL : " << ampl_nlp_object.lagrangian_hessian(iterate) << std::endl;
	std::cout << "Custom : " << problem.lagrangian_hessian(iterate) << std::endl;
		
	return 0;
}