
#include <iostream>
#include <limits>

#include <cmath>

#include "utilities.hh"
#include "line_search.hh"

LineSearchFunction::LineSearchFunction(Nlp &nlp) : FunctionWithNLPState(nlp), penalty_func_(nlp), linear_decrease_func_(nlp) {
	// std::cout << "initializing a linesearcher!" << std::endl;
}
double LineSearchFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) {
	// bugs found here: it_and_step pulled from nlp->initial, which meant it didn't have the right penalty parameter, which caused the linesearch to fail.
	iSQOIterate it_and_step(iterate);
	bool PRINT = false;
	double step_size = 1.0;
	double penfunc_start = penalty_func_(iterate);
	double linear_reduction_start = linear_decrease_func_(iterate, step);
	double machine_precision = std::numeric_limits<double>::epsilon();
	
	double eta = 1e-8;
	double linesearch_step_size_reduction = 0.5;
	if (PRINT) std::cout << std::endl;
	for (size_t alpha_cuts = 0; alpha_cuts < 20; alpha_cuts++) {
		it_and_step.update(iterate, step_size, step);
		double penfunc_step = penalty_func_(it_and_step);
		
		if (PRINT) {
			std::cout << "step_size: " << step_size << std::endl;
			std::cout << "f: " << nlp_->objective(it_and_step) << std::endl;
			std::cout << "c: " << nlp_->constraints_equality(it_and_step) << std::endl;
			std::cout << "alpha_cut: " << alpha_cuts << ": original: " << penfunc_start << "; new: " << penfunc_step << "; reduction: " << (penfunc_start-penfunc_step) << std::endl;
			std::cout << "  rhs: " << - eta*step_size*linear_reduction_start  << " + " << 10*machine_precision*fmax(abs(penfunc_step), 1.0) << " = " << (-eta*step_size*linear_reduction_start+10*machine_precision*std::max(std::abs(penfunc_step), 1.0)) << std::endl;
			std::cout << "mu = " << it_and_step.get_penalty_parameter() << std::endl;
		}
		if (penfunc_step - penfunc_start <= - eta*step_size*linear_reduction_start + 10*machine_precision*fmax(fabs(penfunc_step), 1.0)) {
			break;
		}
		step_size = linesearch_step_size_reduction*step_size;
	}
	return step_size;
}
