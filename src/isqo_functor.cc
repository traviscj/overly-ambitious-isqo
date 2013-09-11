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
#include <iostream>
#include <limits>
#include <string>
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
	
    std::cout << "Hello world! We're serving " << sparse_version_string << " today, please buckle your seat belts and prepare for takeoff..." << std::endl;
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/qpcstair.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
    // DenseAmplNlp problem(problem_file);
    MAGIC_PROBLEM problem(problem_file);
	
	// Utilities for NLP:
	PenaltyFunction penalty_function(problem);
	ResidualFunction kkt_residual(problem);
	SolveQuadraticProgram solve_qp(problem);
	ConstraintViolationFunction constraint_violation(problem);
	LineSearchFunction line_search(problem);
	LinearReductionFunction linear_reduction(problem);
	HessianShifter hessian_shifting_penalty_qp_solve(problem);
	HessianShifter hessian_shifting_feasibility_qp_solve(problem);
	// Other
	TextOutput text_output(problem);
	
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
	bool PRINT=false;
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
	for (iter = 0; iter < maximum_iterations; iter ++ ) {
        // std::cout << std::endl << "penalty: " << penalty_iterate.penalty_parameter_ << std::endl;
		iSQOStep combination_step(problem.num_primal(), problem.num_dual_eq(), problem.num_dual_ieq(), -42, -42);
		
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
		
		//////////////////////////
		// ALGORITHM A // STEP 2
		//////////////////////////
		if ((kkt_residual(penalty_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) < termination_threshold)) {
			std::cout << std::endl << "Termination 2a - optimality" << std::endl;
			break;
		}
		if ((kkt_residual(feasibility_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) > termination_threshold) && (penalty_iterate.get_penalty_parameter() < termination_penalty_threshold)) {
			std::cout << std::endl << "Termination 2b - infeasible problem!" << std::endl;
			break;
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 3
		//////////////////////////
		// Penalty problem is set up AND SOLVED, 
		MAGIC_SUBPROBLEM penalty_subproblem(problem, penalty_iterate);
		iSQOStep penalty_step = hessian_shifting_penalty_qp_solve(penalty_subproblem);
		iSQOStep feasibility_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq(), -43, -43);
		
		// per-iteration variable setup.
		std::string steptype = "4a";
		double combination_step_contribution_from_penalty_step = -1.0;
		double next_penalty_parameter;
		if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*constraint_violation(penalty_iterate) + 10*machine_precision*linear_reduction(penalty_iterate,penalty_step)) {
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
			MAGIC_SUBPROBLEM feasibility_subproblem(problem, feasibility_iterate);
			feasibility_step = hessian_shifting_feasibility_qp_solve(feasibility_subproblem);
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
                    // assert(! (num_comb_reductions >= max_num_comb_reductions));
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
                    assert(potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator > 0);
					next_penalty_parameter = std::min(  penalty_parameter_reduction_factor*penalty_iterate.get_penalty_parameter(), 
                                                                    potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator);
				}
			}

			// combination gets the convex combination of penalty and feasibility:
			combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
			
			text_output.subproblem(hessian_shifting_feasibility_qp_solve.get_last_shift(), feasibility_iterate, feasibility_subproblem, feasibility_step);
		}
		text_output.subproblem(hessian_shifting_penalty_qp_solve.get_last_shift(), penalty_iterate, penalty_subproblem, penalty_step);
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, combination_step_contribution_from_penalty_step);
        
        penalty_iterate.set_penalty_parameter(next_penalty_parameter);
		double step_size = line_search(penalty_iterate, combination_step);
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
	}
	
	if (iter == maximum_iterations) {
		std::cout << "Failure - Did not converge" << std::endl;
	}
	std::cout << penalty_iterate.get_primal_values() << std::endl;
	return 0;
}


