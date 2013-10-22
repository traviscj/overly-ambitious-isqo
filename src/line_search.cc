// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#include <iostream>
#include <limits>

#include <cmath>

#include "isqo_config.hh"
#include "utilities.hh"
#include "line_search.hh"

LineSearchFunction::LineSearchFunction(iSQOControlPanel &control, Nlp &nlp) : FunctionWithNLPState(control, nlp), penalty_func_(control, nlp), linear_decrease_func_(control, nlp) {
	// std::cout << "initializing a linesearcher!" << std::endl;
}
double LineSearchFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) {
	// bugs found here: it_and_step pulled from nlp->initial, which meant it didn't have the right penalty parameter, which caused the linesearch to fail.
    std::cout << "**** Line Search init." << std::endl;
	iSQOIterate it_and_step(iterate);
	bool PRINT = false;
	double step_size = 1.0;
	double penfunc_start = penalty_func_(iterate);
    // THIS WAS THE ORIGINAL LINE:
    // double linear_reduction_start = linear_decrease_func_(iterate, step);
    double linear_reduction_start = linear_decrease_func_(iterate, step);
	double machine_precision = std::numeric_limits<double>::epsilon();
	
	double eta = 1e-8;
	double linesearch_step_size_reduction = 0.5;
    std::cout << "**** Starting Line Search." << std::endl;
    
    size_t alpha_cuts = 0;
	for (; alpha_cuts < 20; alpha_cuts++) {
		it_and_step.update(iterate, step_size, step);
		double penfunc_step = penalty_func_(it_and_step);
		
		if (PRINT) {
			std::cout << "step_size: " << step_size << std::endl;
			std::cout << "f: " << nlp_->objective(it_and_step) << std::endl;
			std::cout << "c: " << nlp_->constraints_equality(it_and_step) << std::endl;
			std::cout << "alpha_cut: " << alpha_cuts << ": original: " << penfunc_start << "; new: " << penfunc_step << "; reduction: " << (penfunc_start-penfunc_step) << std::endl;
			std::cout << "  rhs: " << - eta*step_size*linear_reduction_start  << " + " << 10*machine_precision*fmax(fabs(penfunc_step), 1.0) << " = " << (-eta*step_size*linear_reduction_start+10*machine_precision*fmax(fabs(penfunc_step), 1.0)) << std::endl;
			std::cout << "mu = " << it_and_step.get_penalty_parameter() << std::endl;
		}
        
        std::cout << "***** alpha: " << step_size << "; pf step: " << penfunc_step << "; "
                    << "pf start: " << penfunc_start << "; "
                    << "lredstart: " << linear_reduction_start << "; "
                    << "(entire) relative part: " << 10*machine_precision*fmax(fabs(penfunc_step), 1.0)
                    << std::endl;
		if (penfunc_step - penfunc_start <= - eta*step_size*linear_reduction_start + 10*machine_precision*fmax(fabs(penfunc_step), 1.0)) {
			break;
		}
		step_size = linesearch_step_size_reduction*step_size;
	}
    std::cout << "**** Finished Line Search; alpha: " << step_size << "; alpha_cuts: " << alpha_cuts << std::endl;
	return step_size;
}
