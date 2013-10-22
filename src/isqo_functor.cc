// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
// overly-ambitious-iSQO
//  - Travis C. Johnson (traviscj@traviscj.com)
//  - August 2013

// basic philosophy is:
// - keep main loop as clean and as close to algo in the paper as humanly possible, at the expense of anything else.
// - don't try to fit functions (e.g. linear reduction) into classes for iterates
//   - INSTEAD: write functions that stand on their own (but functors, so they can have some persistent state)
// - don't try to fit multiple types of iterates into one class
//   - INSTEAD: separate the iterates from steps from subproblems from hessian shifting from...
// - don't try to fit more functionality than absolutely necessary into classes.
//   - INSTEAD: have easy-to-debug classes with a single simple point: handling hessian shifts, for instance.
// - don't try to pass in the kitchen sink simply to be able to evaluate some basic things.
//   - INSTEAD: Pass the NLP class, which returns f(x), c(x), \bar{c}(x), take it from there.
// - DO NOT use short (1-2 character) or greek symbols as variable names.
//   - INSTEAD: spell out exactly what you are iterating over. Primal indices like "i"? Try primal_index. Linesearch step size called "alpha"? Try linesearch_step_size.
//   - Exception: External class interfaces are a bit messier.
//   - Gotcha: A bit more consistent naming scheme would improve things even more.
// - whenever possible, return ONE object by value. 
//   - INSTEAD: let C++11x worry about RVO. 
//   - Exception: iSQOStep class... but maybe it doesn't need to be?
// - whenever possible, pass required objects by const reference.
// - whenever possible, use std::vector<double> instead of any other std::vectors.
// - for unfinished features, use assert to ensure they are not used.

// future feature list:
// - DONE !!! variable bounds !!!
// - DONE !!! more robust shifting strategy and/or Quasi-Newton. (Quasi-Newton might still be an interesting experiment)
// - DONE !!! rewrite some of the matrix computations into the matrix class. ) !!!
// - DONE !!! integrate DJ/CK/AW/et al's work on sparse qpOASES. !!!
// - constraint/objective scaling (should fix HS99... )
// - other iSQO features (algo II-features.)
// - iQP interface! -- but how to do fallback/etc?
// - logbook-style reporting/logging, with sql.
// - use qpoases matrices instead of custom matrix code
// - sparse matrices from ampl/to qpOASES.
// - BQPD interface?
// - read parameters from files, store in some structure.
// - second order correction (might fix HS067...)
// - use multiple qp objects for getting subproblem solutions.
// - const all applicable functions.
// - stay bound feasible option? (eg, dallas* need feasible primal variables.)
// - 

// justification for weirder choices:
// - I use "functors" for some of the paper-defined functions, to capture program state. They work by overloading operator(), which is a big weird, but they are essentially functions with state. One key thing: State is still minimized in them--most are just wrappers around basic NLP functions. With this I was trying to accomplish clear input-output relations, instead of inventing functions that do incremental changes to a huge block of state. One downside to this, as written, is that many things get re-evaluated(This could be minimized by some caching of )

// future ideas:
// - would like to extend the FunctionWithNLPState class to cache the contents of the arguments, and return the last value if they haven't changed. I think this could be somewhat straightforward, but maybe I'm underestimating. This is another major "point" to doing 
// 

// helpful hints:
// - OSX: for alloc bugs, run with: export DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib. On other systems, valgrind.
// - currently, here, the reported KKT is the two-norm... in the matlab implementation, it is the inf norm, scaled by a bunch of things.

// remaining potential bugs:
// hs067.out:Failure - Did not converge - tiny step length... maybe bad step from QP solver? -- might be best fixed by a second-order correction.
// hs099.out:Failure - Did not converge - this *should* be fixed with problem rescaling - should add problem rescaling.
// other HS failures:
// hs087.out:Failure - Did not converge - was also iter limit in MATLAB implmentation.
// hs098.out:Failure - Did not converge - subproblem failure in MATLAB implementation.
// hs106.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs109.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs114.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs117.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs99exp.out:Failure - Did not converge - iteration limit in MATLAB implementation.

// standard libraries:
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>

// iSQO headers:
#include "isqo.hh"

#ifdef ISQO_USE_DENSE
const std::string sparse_version_string = "DenseAmplNlp and iSQOQuadraticSubproblem";
#define MAGIC_PROBLEM DenseAmplNlp
#define MAGIC_SUBPROBLEM iSQOQuadraticSubproblem
#else
const std::string sparse_version_string = "SparseAmplNlp and iSQOSparseQuadraticSubproblem";
#define MAGIC_PROBLEM SparseAmplNlp
#define MAGIC_SUBPROBLEM iSQOQuadraticSubproblem
#endif

int main(int argc, char **argv) {
    
    iSQOControlPanel control("config_reader_input");
    control.update_settings();
    control.print();
    
    std::cout << "Hello world! We're serving " << sparse_version_string << " today, please buckle your seat belts and prepare for takeoff..." << std::endl;
    std::cout << "* AMPL Setup Caterwalling..." << std::endl;
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/qpcstair.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
    // DenseAmplNlp problem(problem_file);
    MAGIC_PROBLEM problem(problem_file);
	
	// Utilities for NLP:
	PenaltyFunction penalty_function(control, problem);
	ResidualFunction kkt_residual(control, problem);
    // SolveQuadraticProgram solve_qp(control, problem);
	ConstraintViolationFunction constraint_violation(control, problem);
	LineSearchFunction line_search(control, problem);
	LinearReductionFunction linear_reduction(control, problem);
	HessianShifter hessian_shifting_penalty_qp_solve(control, problem);
	HessianShifter hessian_shifting_feasibility_qp_solve(control, problem);
	// Other
	TextOutput text_output(control, problem);
	
    text_output.nlp();
    
	// Paper-defined Parameter Values
	double linear_decrease_threshold = 1e-1; // epsilon in the paper.
	double linear_reduction_threshold_for_penalty_reduction = 1e-2; // beta in the paper.
	double mostly_feasibility_step_threshold = 1e-3; // tau in the paper.
	double penalty_parameter_reduction_factor = 0.2; // delta in paper.
	double termination_threshold = 1e-6;
	double termination_penalty_threshold = 1e-8;
	double hessian_min_convexity = 1e-8;	// theta in the paper.
	// Code-defined parameter values:
	double convex_combination_search_decrease = .9;
    // bool PRINT=false;
	size_t maximum_iterations = 1000;
    size_t max_num_comb_reductions = 100;
	double machine_precision = std::numeric_limits<double>::epsilon();
	
	text_output.start();
	
	//////////////////////////
	// ALGORITHM A // STEP 1
	//////////////////////////
	iSQOIterate penalty_iterate = problem.initial(1.0e-1);
	iSQOIterate feasibility_iterate = problem.initial(0.0);

	size_t iter=-1;
    std::cout << "* Main iteration Loop: " << std::endl;
    
	for (iter = 0; iter < maximum_iterations; iter ++ ) {
        if (iter == 0)
            std::cout << "** front matter..." << std::endl;
        else
            std::cout << "*** front matter..." << std::endl;
        // std::cout << std::endl << "penalty: " << penalty_iterate.penalty_parameter_ << std::endl;
        double pre_objective = problem.objective(penalty_iterate);
        double pre_infeas = constraint_violation(penalty_iterate);
        double pre_merit = penalty_function(penalty_iterate);
        double pre_feaskkt = kkt_residual(feasibility_iterate);
        double pre_penkkt = kkt_residual(penalty_iterate);
		iSQOStep combination_step(problem.num_primal(), problem.num_dual_eq(), problem.num_dual_ieq(), -42, -42);
		
        std::cout << "** Start of iteration #" << iter << "; " ;
            // << std::endl;
        std::cout << std::scientific << std::setprecision(2) << std::showpos;
        std::cout << "obj: " << pre_objective << "; ";
        std::cout << "infeas: " << pre_infeas << "; ";
		std::cout << "pen: " << penalty_iterate.get_penalty_parameter() << "; ";
        std::cout << "merit: " << pre_feaskkt << "; ";
        std::cout << "feaskkt: " << pre_penkkt << "; ";
        std::cout << "penkkt: " << pre_penkkt << "; ";
        std::cout << std::endl;
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
        
		//////////////////////////
		// ALGORITHM A // STEP 2
		//////////////////////////
		if ((kkt_residual(penalty_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) < termination_threshold)) {
			std::cout << std::endl << "* Final Status: Termination 2a - optimality" << std::endl;
            text_output.finish_success_opt();
			break;
		}
		if ((kkt_residual(feasibility_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) > termination_threshold) && (penalty_iterate.get_penalty_parameter() < termination_penalty_threshold)) {
			std::cout << std::endl << "* Final Status: Termination 2b - infeasible problem!" << std::endl;
            text_output.finish_success_inf();
			break;
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 3
		//////////////////////////
		// Penalty problem is set up AND SOLVED, 
        std::cout << "*** creating the PENALTY subproblem: "<< std::endl;
        std::cout << "**** penalty_iterate: " << penalty_iterate << std::endl;
		MAGIC_SUBPROBLEM penalty_subproblem(control, problem, penalty_iterate);
        // std::cout << "Penalty QP start:" << std::endl;
		iSQOStep penalty_step = hessian_shifting_penalty_qp_solve(penalty_iterate, penalty_subproblem);
        // std::cout << "Penalty QP end:" << std::endl;
		iSQOStep feasibility_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq(), -43, -43);
		
		// per-iteration variable setup.
		std::string steptype = "4a";
		double combination_step_contribution_from_penalty_step = -1.0;
		double next_penalty_parameter;

        std::cout << "*** matrices for shift: " << hessian_shifting_penalty_qp_solve.get_last_shift() << std::endl;
        std::shared_ptr<matrix_base_class> nlp_hessian(problem.lagrangian_hessian(penalty_iterate));
        
        std::cout << "**** NLP hessian: " << std::endl << nlp_hessian << std::endl;
        std::cout << "**** subproblem hessian: " << std::endl << penalty_subproblem.hessian_ << std::endl;
        std::cout << "*** Output for subsequent KKT residual call:" << std::endl;
        double pen_kkt_resid = kkt_residual(penalty_iterate, penalty_subproblem, penalty_step);
        std::cout << "*** kkt_residual(penalty_iterate,penalty_step): " << pen_kkt_resid << std::endl;
        
        double val = linear_reduction(penalty_iterate,penalty_step);
        std::cout << "*** linear_reduction(penalty_iterate,penalty_step): " << val << std::endl;
        assert(val>=-1e-12); // need val to be... mostly positive.
		if (val >= linear_decrease_threshold*constraint_violation(penalty_iterate) + 10*machine_precision*val) {
			//////////////////////////
			// ALGORITHM A // STEP 3a
			combination_step_contribution_from_penalty_step = 1.0;
			text_output.subproblem_skip();
            
            next_penalty_parameter = penalty_iterate.get_penalty_parameter();
			combination_step.set_primal(penalty_step);
		} else {
			
			//////////////////////////
			// ALGORITHM A // STEP 4
			//////////////////////////
			// Feasibility problem is set up and solved.
            // std::cout << "Feasibility QP start:" << std::endl;
            std::cout << "*** creating the FEASIBILITY subproblem: "<< std::endl;
			MAGIC_SUBPROBLEM feasibility_subproblem(control, problem, feasibility_iterate);
			feasibility_step = hessian_shifting_feasibility_qp_solve(feasibility_iterate, feasibility_subproblem);
            std::cout << "*** linear_reduction(penalty_iterate,feasibility_step): " << linear_reduction(penalty_iterate,feasibility_step) << std::endl;
            
            // std::cout << "Feasibility QP end:" << std::endl;
            assert(feasibility_step.get_status() == 0);
			
			if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*linear_reduction(feasibility_iterate,feasibility_step)) {
				//////////////////////////
				// ALGORITHM A // STEP 4a
				combination_step_contribution_from_penalty_step = 1.0;
				steptype = "5a";
			} else {
				//////////////////////////
				// ALGORITHM A // STEP 4b
				combination_step_contribution_from_penalty_step = 1.0;
				
				int num_comb_reductions = 0;
				combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
				while (! (linear_reduction(feasibility_iterate,combination_step)  >=linear_decrease_threshold*linear_reduction(feasibility_iterate,feasibility_step) + 10*machine_precision*linear_reduction(feasibility_iterate,combination_step)) ) {
					combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
					
					combination_step_contribution_from_penalty_step = convex_combination_search_decrease*combination_step_contribution_from_penalty_step;
					++num_comb_reductions;
                    if (num_comb_reductions >= max_num_comb_reductions) {
                        combination_step_contribution_from_penalty_step = 0;
                        combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
                        break;
                    }
                    // assert(! );
				}
				steptype = "5b";
				
				// beta: 
				double linear_decrease_in_penalty_combination = linear_reduction(penalty_iterate,combination_step);
				double linear_decrease_in_feasibility_combination = linear_reduction(feasibility_iterate,combination_step);
				if (linear_decrease_in_penalty_combination >= linear_reduction_threshold_for_penalty_reduction*linear_decrease_in_feasibility_combination + 10*machine_precision*linear_decrease_in_penalty_combination) {
					if (combination_step_contribution_from_penalty_step >= mostly_feasibility_step_threshold) {
						// no-op: keep the current penalty parameter.
                        next_penalty_parameter = penalty_iterate.get_penalty_parameter();
					} else {
                        // penalty_iterate.set_penalty_parameter();
                        next_penalty_parameter = penalty_parameter_reduction_factor*penalty_iterate.get_penalty_parameter();
					}
				} else {
					std::vector<double> gradient = problem.objective_gradient(penalty_iterate);
					double potential_new_penalty_parameter_numerator = (1-linear_reduction_threshold_for_penalty_reduction)*linear_decrease_in_feasibility_combination;
					double potential_new_penalty_parameter_denominator = combination_step.x_dot_product(gradient) + hessian_min_convexity*pow(combination_step.x_norm(),2);
                    // assert(potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator > 0);
					next_penalty_parameter = std::min(  penalty_parameter_reduction_factor*penalty_iterate.get_penalty_parameter(), 
                                                                    INFINITY*potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator);
				}
			}

			// combination gets the convex combination of penalty and feasibility:
			combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
			
			text_output.subproblem(hessian_shifting_feasibility_qp_solve.get_last_shift(), feasibility_iterate, feasibility_subproblem, feasibility_step);
		}
		text_output.subproblem(hessian_shifting_penalty_qp_solve.get_last_shift(), penalty_iterate, penalty_subproblem, penalty_step);
        std::cout << "*** Termination Type: " << steptype << std::endl;
        
        double val2 = linear_reduction(penalty_iterate,combination_step);
        std::cout << "*** linear_reduction(penalty_iterate,combination_step): " << val2 << std::endl;
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
        
        std::cout << "*** pre .post " << std::endl;
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, combination_step_contribution_from_penalty_step);
        std::cout << "*** post .post " << std::endl;
        
        penalty_iterate.set_penalty_parameter(next_penalty_parameter);
        std::cout << "*** pre .line_search " << std::endl;
		double step_size = line_search(penalty_iterate, combination_step);
        std::cout << "*** post .line_search; step size: " << step_size << std::endl;
        text_output.line_search(step_size);
		
		//////////////////////////
		// ALGORITHM A // STEP 6
		//////////////////////////
		// update penalty iterate & dual values:
		// update FEASIBILITY iterate first, since penalty iterate update changes this answer.
		penalty_iterate.update(penalty_iterate, step_size, combination_step);
		penalty_iterate.update_dual(penalty_step);

        // feasibility_iterate.update(penalty_iterate, step_size, combination_step);
        feasibility_iterate.assign_primal(*(penalty_iterate.get_primal_values()));
		feasibility_iterate.update_dual(feasibility_step);
        
        problem.reset_cache();
        hessian_shifting_penalty_qp_solve.save_qp_state();
        if (feasibility_step.get_status() == 0) hessian_shifting_feasibility_qp_solve.save_qp_state();
        
        std::cout << hessian_shifting_penalty_qp_solve.get_info_str();
        std::cout << hessian_shifting_feasibility_qp_solve.get_info_str();
        // std::cout << combination_info.str();
	}
	
	if (iter == maximum_iterations) {
        std::cout << "* Final Status: Failure - Did not converge" << std::endl;
        text_output.finish_fail();
	}
    // std::cout << penalty_iterate.get_primal_values() << std::endl;
	return 0;
}


