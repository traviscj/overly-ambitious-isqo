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
// - DONE... rewrite some of the matrix computations into the matrix class. (But still need to check for stragglers...)
// - constraint/objective scaling (should fix HS99... )
// - other iSQO features (algo II-features.)
// - iQP interface! -- but how to do fallback/etc?
// - integrate DJ/CK/AW/et al's work on sparse qpOASES.
// - logbook-style reporting/logging, with sql.
// - use qpoases matrices instead of custom matrix code
// - sparse matrices from ampl/to qpOASES.
// - BQPD interface?
// - read parameters from files, store in some structure.
// - second order correction (might fix HS067...)
// - use multiple qp objects for getting subproblem solutions.
// - const all applicable functions.
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
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

// 
#include <qpOASES.hpp>

// iSQO headers:
#include "utilities.hh"
#include "step.hh"
#include "iterate.hh"
#include "nlp.hh"
#include "nlp_hs014.hh"
#include "nlp_ampl.hh"
#include "nlp_state.hh"
#include "constraint_violation.hh"
#include "penalty_function.hh"
#include "subproblem.hh"
#include "residual_function.hh"
#include "solve_subproblem.hh"
#include "linear_model.hh"
#include "linear_model_reduction.hh"
#include "line_search.hh"
#include "hessian_shifter.hh"

class TextOutput : public FunctionWithNLPState {
public:
	TextOutput (Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp),linear_decrease_func_(nlp), pen_func_(nlp), residual_func_(nlp) {
		
	}
	
	void start() {
		printf("-----|");
		printf("----------------------|");
		printf("----------------------|");
		printf("----------------------&");
		printf("-----------------[ FEASIBILITY ]-----------------|");
		printf("-------------------[ PENALTY ]-------------------|");
		printf("--------------------------------------------|");
		printf("----------\n");
		printf("%s",output_desc_pre_);
		printf("%s",output_desc_subprob_);
		printf("%s",output_desc_subprob_);
		printf("%s",output_desc_post_);
		printf("-----|");
		printf("----------------------|");
		printf("----------------------|");
		printf("----------------------&");
		printf("-------------------------------------------------|");
		printf("-------------------------------------------------|");
		printf("--------------------------------------------|");
		printf("----------\n");
	}
	void pre(size_t iter, const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate) const {
		// iter, problem, feas iterate, pen iterate, constraintviolation, iterate, penfunc, residual_func
		//  -> params: iter, feas iter, pen iter
		//	-> state: nlp_, constraint_violation_func_, pen_func_, residual_func_
		printf(output_format_pre_, 
				iter, 
				nlp_->objective(penalty_iterate), constraint_violation_func_(penalty_iterate),
				penalty_iterate.penalty_parameter_, pen_func_(penalty_iterate),
				residual_func_(feasibility_iterate), residual_func_(penalty_iterate)
				);
	}
	void subproblem(double shift, const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// shift, step, linear_decrease_func, residual_func
		//  -> params: shift, step
		//	-> state: linear_decrease_func, residual_func_
		printf(output_format_subprob_,
				shift, step.status_, step.x_norm(),
				linear_decrease_func_(iterate, step), residual_func_(iterate, subproblem, step)
					);
	}
	void subproblem_skip() {
		printf("         -   -           -          -          - |");
	}
	void post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, std::string step_type, double step_mix, double alpha) {
		// step, feas_iter, penalty_iter, alpha
		//  -> params: step, feas_iter, penalty_iter, alpha
		printf(output_format_post_,
				step_type.c_str(),
				step_mix, combination_step.x_norm(), linear_decrease_func_(feasibility_iterate,combination_step), linear_decrease_func_(penalty_iterate, combination_step),
				alpha
					);
	}
protected:
	// Nlp *nlp_;
	ConstraintViolationFunction constraint_violation_func_;
	LinearReductionFunction linear_decrease_func_;
	PenaltyFunction pen_func_;
	ResidualFunction residual_func_;
	
	static const char output_desc_pre_[];
	static const char output_desc_subprob_[];
	static const char output_desc_post_[];
	static const char output_format_pre_[];
	static const char output_format_subprob_[];
	static const char output_format_post_[];
};
const char TextOutput::output_desc_pre_[] = " it  |       obj     infeas |       pen      merit |   feaskkt     penkkt &";
const char TextOutput::output_desc_subprob_[] = "     shift  msg     ||d||     penred        res  |";
const char TextOutput::output_desc_post_[] = " TT   CvxComb    ||d||   FeasRed     PenRed |    alpha\n";
const char TextOutput::output_format_pre_[] = " %3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
const char TextOutput::output_format_subprob_[] = " %9.2e  %3d  %9.2e  %9.2e  %9.2e |";
const char TextOutput::output_format_post_[] = " %s %9.2e %9.2e %9.2e %9.2e |%9.2e\n";


int main(int argc, char **argv) {
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs067.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
	AmplNlp problem(problem_file);
	iSQOIterate penalty_iterate = problem.initial();
    // penalty_iterate.primal_values_[0] = 1.25*3;
    // penalty_iterate.primal_values_[1] = 1.25*3;
    // penalty_iterate.primal_values_[2] = 1.25*3;
    // penalty_iterate.primal_values_[3] = 1.25*3;
    // penalty_iterate.primal_values_[4] = 1.25*3;
    // penalty_iterate.primal_values_[5] = 1.25*3;
    // penalty_iterate.primal_values_[7] = 1.25*3;
    // penalty_iterate.primal_values_[8] = 1.25*3;
    // penalty_iterate.primal_values_[9] = 1.25*3;
    
	
    std::cout << "objective: " << problem.objective(penalty_iterate) << std::endl;
	std::cout << "constraints: "<< problem.constraints_equality(penalty_iterate) << std::endl;
	
    qpOASES::SparseMatrix sparse_jacobian = problem.constraints_equality_jacobian_sparse(penalty_iterate);
    // std::cout << "eq jacobian: " << ) << std::endl;
    
}
int main3() {
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
int main2() {
	std::cout << "okay, so..." << std::endl;
	
	// tests for constraint conversation:
	// equality constrained problems: (to c(x) = 0)
	//	* 0 <= c_1(x) <= 0  ---> A_1d = 0 - c_1(x_k)
	//	* 1 <= c_2(x) <= 1  ---> A_2d = 1 - c_2(x_k)
	// inequality constrained problems: (to barc(x) <= 0)
	//  * double constrained: 	-1 <= c_3(x) <= 5  --> A_3d <= 5 - c_3(x_k) AND -A_3d <= -(-1) - (-c_3(x_k))
	//			c_3(x) <= 5  ==> 
	//			-1 <= c_3(x) ==> -c_3(x) <= -(-1) ==> -A_3d - c_3(x_k) <= -(-1) ==> -A_3d  <= -(-1) - (-c_3(x_k))
	//	* double constrained: 	0 <= c_4(x) <= 5 --> A_4d <= 5 - c_4(x_k) AND -A_4d <= -(0) - (-c_4(x_k))
	//	* double constrained: 	-1 <= c_5(x) <= 0 --> A_5d <= 0 - c_5(x_k) AND -A_5d <= -(-1) - (-c_5(x_k))
	//	* only upper: -INF <= c_6(x) <= 0 --> A_6d <= 0 - c_6(x_k)
	//	* only upper: -INF <= c_7(x) <= 5 --> A_7d <= 5 - c_7(x_k)
	//	* only lower: 0 <= c_8(x) <= INF --> -A_8d <= -(0) - c_8(x_k)
	//	* only lower: 5 <= c_9(x) <= INF --> -A_9d <= -(5) - c_9(x_k)
	
	// 
	// tests for variable conversion:
	// 
	// AmplNlp ampl_nlp_object("/Users/traviscj/optimization/cute_nl_nopresolve/hs016.nl");
	AmplNlp ampl_nlp_object("/Users/traviscj/ampl_bound_reformulation.nl");
	iSQOIterate iterate = ampl_nlp_object.initial();
	std::cout << "initial iterate: " << iterate.primal_values_ << std::endl;
	// std::cout << "Zeroth Order:" << std::endl;
	std::cout << "===[Objective]==============" << std::endl;
	double ampl_obj = ampl_nlp_object.objective(iterate);
	assert(ampl_obj == 5.921590000000000e+00);//, "ampl_obj wrong");
	std::cout << "  AMPL : " << ampl_nlp_object.objective(iterate) << std::endl;
	
	std::cout << "===[Objective Gradient]==============" << std::endl;
	std::vector<double> ampl_obj_gradient = ampl_nlp_object.objective_gradient(iterate);
	assert(ampl_obj_gradient[0] == 1.0e0);
	assert(ampl_obj_gradient[1] == 1.0e0);
	std::cout << "  AMPL : " << ampl_obj_gradient << std::endl;
	
	
	std::cout << "===[Equality Constraints]==============" << std::endl;
	std::vector<double> con_eq_values = ampl_nlp_object.constraints_equality(iterate);
	std::cout << "  AMPL : " << con_eq_values << std::endl;
	assert(con_eq_values[0]==-41.924375456199996); // negative because switches var to left-hand-side.
	assert(con_eq_values[1]==141.09951409669998);
	
	std::cout << "===[Equality Jacobian]==============" << std::endl;
	matrix eq_jacobian = ampl_nlp_object.constraints_equality_jacobian(iterate);
	std::cout << "  AMPL : " << eq_jacobian << std::endl;
	assert(eq_jacobian.get(0,0) == -1.256636000000000e+01);
	assert(eq_jacobian.get(0,1) == -1.668000000000000e+01);
	assert(eq_jacobian.get(1,0) == 4.398226000000000e+01);
	assert(eq_jacobian.get(1,1) == 6.116000000000000e+01);
	
	
	
	
	std::cout << "===[Inequality Constraints]==============" << std::endl;
	std::vector<double> ieq_values = ampl_nlp_object.constraints_inequality(iterate);
	std::cout << "  AMPL : " << ieq_values << std::endl;
	double ieq_tol = 1e-12;
	assert(abs(ieq_values[0]) <= (abs(-3.482753668339000e+02) + ieq_tol*abs(ieq_values[0])));
	assert(abs(ieq_values[1]) <= (abs(-6.820391459396999e+02) + ieq_tol*abs(ieq_values[1])));
	assert(abs(ieq_values[2]) <= (abs(3.362753668339000e+02) + ieq_tol*abs(ieq_values[2])));
	assert(abs(ieq_values[3]) <= (abs(7.346270723082998e+02) + ieq_tol*abs(ieq_values[3])));
	
	// -6.820391459396999e+02
	
	std::cout << "===[Inequality Jacobian]==============" << std::endl;
	matrix ieq_jacobian = ampl_nlp_object.constraints_inequality_jacobian(iterate);
	std::cout << "  AMPL : " << ieq_jacobian << std::endl;
	assert_close(ieq_jacobian.get(0,0), -1.193804200000000e+02, 1e-10); // c2 lower
	assert_close(ieq_jacobian.get(0,1), -1.278800000000000e+02, 1e-10);
	assert_close(ieq_jacobian.get(1,0), -2.324776600000000e+02, 1e-10); // c3 lower
	assert_close(ieq_jacobian.get(1,1), -2.279600000000000e+02, 1e-10);
	assert_close(ieq_jacobian.get(2,0), 1.193804200000000e+02, 1e-10);  // c2 upper
	assert_close(ieq_jacobian.get(2,1), 1.278800000000000e+02, 1e-10);  
	assert_close(ieq_jacobian.get(3,0), 2.701767400000000e+02, 1e-10);  // c4 upper
	assert_close(ieq_jacobian.get(3,1), 2.613200000000000e+02, 1e-10);
	
	
	
	
	
	
	// assert(ieq_jacobian.get(3,0) == 1.193804200000000e+02); assert(ieq_jacobian.get(3,1) == 1.278800000000000e+02);  // c2 upper
	
	std::cout << std::endl;
	std::cout << "TESTED FRONTIER!" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	
	std::cout << "===[Lagrangian Hessian]==============" << std::endl;
	matrix lh_f_only = ampl_nlp_object.lagrangian_hessian(iterate);
	std::cout << "  AMPL : " << lh_f_only << std::endl;
	assert_close(lh_f_only.get(0,0), 0.00e0, 1e-10);  // c2 upper
	assert_close(lh_f_only.get(0,1), 0.00e0, 1e-10);  
	assert_close(lh_f_only.get(1,0), 0.00e0, 1e-10);  // c4 upper
	assert_close(lh_f_only.get(1,1), 0.00e0, 1e-10);
	
	iterate.dual_eq_values_[0] = 1.0;
	matrix lh_c0_only = ampl_nlp_object.lagrangian_hessian(iterate);
	std::cout << "  AMPL : " << lh_c0_only << std::endl;
	iterate.dual_eq_values_[0] = -1.0;
	matrix lh_c0_only_flip = ampl_nlp_object.lagrangian_hessian(iterate);
	std::cout << "  AMPL : " << lh_c0_only_flip << std::endl;
	
	
	return 0;
}
int main1() {
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
int main0(int argc, char **argv) {
	
	std::string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs067.nl");
	if (argc>1) {
		problem_file = std::string(argv[1]);
	}
	AmplNlp problem(problem_file);
	
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
	iSQOIterate penalty_iterate = problem.initial();
	penalty_iterate.penalty_parameter_ = 1.0e-1;

	iSQOIterate feasibility_iterate = problem.initial();
	feasibility_iterate.penalty_parameter_ = 0.0;
	size_t iter=-1;
	for (iter = 0; iter < maximum_iterations; iter ++ ) {
		iSQOStep combination_step(problem.num_primal(), problem.num_dual_eq(), problem.num_dual_ieq());
		
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
		
		//////////////////////////
		// ALGORITHM A // STEP 2
		//////////////////////////
		if ((kkt_residual(penalty_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) < termination_threshold)) {
			std::cout << std::endl << "Termination 2a - optimality" << std::endl;
			break;
		}
		if ((kkt_residual(feasibility_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) > termination_threshold) && (penalty_iterate.penalty_parameter_ < termination_penalty_threshold)) {
			std::cout << std::endl << "Termination 2b - infeasible problem!" << std::endl;
			break;
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 3
		//////////////////////////
		// Penalty problem is set up AND SOLVED, 
		iSQOQuadraticSubproblem penalty_subproblem(problem, penalty_iterate);
		iSQOStep penalty_step = hessian_shifting_penalty_qp_solve(penalty_iterate, penalty_subproblem);
		iSQOStep feasibility_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());
		
		// per-iteration variable setup.
		std::string steptype = "4a";
		double combination_step_contribution_from_penalty_step = -1.0;
		
		if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*constraint_violation(penalty_iterate) + 10*machine_precision*linear_reduction(penalty_iterate,penalty_step)) {
			//////////////////////////
			// ALGORITHM A // STEP 3a
			combination_step_contribution_from_penalty_step = 1.0;
			text_output.subproblem_skip();

			combination_step.set_primal(penalty_step);
		} else {
			
			//////////////////////////
			// ALGORITHM A // STEP 4
			//////////////////////////
			// Feasibility problem is set up and solved.
			iSQOQuadraticSubproblem feasibility_subproblem(problem, feasibility_iterate);
			feasibility_step = hessian_shifting_feasibility_qp_solve(feasibility_iterate, feasibility_subproblem);
			
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
					assert(! (num_comb_reductions >= max_num_comb_reductions));
				}
				steptype = "5b";
				
				// beta: 
				double linear_decrease_in_penalty_combination = linear_reduction(penalty_iterate,combination_step);
				double linear_decrease_in_feasibility_combination = linear_reduction(feasibility_iterate,combination_step);
				if (linear_decrease_in_penalty_combination >= linear_reduction_threshold_for_penalty_reduction*linear_decrease_in_feasibility_combination + 10*machine_precision*linear_decrease_in_penalty_combination) {
					if (combination_step_contribution_from_penalty_step >= mostly_feasibility_step_threshold) {
						// no-op: keep the current penalty parameter.
					} else {
						penalty_iterate.penalty_parameter_ = penalty_parameter_reduction_factor*penalty_iterate.penalty_parameter_;
					}
				} else {
					std::vector<double> gradient = problem.objective_gradient(penalty_iterate);
					double potential_new_penalty_parameter_numerator = (1-linear_reduction_threshold_for_penalty_reduction)*linear_decrease_in_feasibility_combination;
					double potential_new_penalty_parameter_denominator = combination_step.x_dot_product(gradient) + hessian_min_convexity*pow(combination_step.x_norm(),2);
					penalty_iterate.penalty_parameter_ = std::min(penalty_parameter_reduction_factor*penalty_iterate.penalty_parameter_, potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator);
				}
			}

			// combination gets the convex combination of penalty and feasibility:
			combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
			
			text_output.subproblem(hessian_shifting_feasibility_qp_solve.get_last_shift(), feasibility_iterate, feasibility_subproblem, feasibility_step);
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
		double step_size = line_search(penalty_iterate, combination_step);
		
		text_output.subproblem(hessian_shifting_penalty_qp_solve.get_last_shift(), penalty_iterate, penalty_subproblem, penalty_step);
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, combination_step_contribution_from_penalty_step, step_size);
		
		//////////////////////////
		// ALGORITHM A // STEP 6
		//////////////////////////
		// update penalty iterate & dual values:
		penalty_iterate.update(penalty_iterate, step_size, combination_step);
		penalty_iterate.update_dual(penalty_step);
		// update feasibility iterate & dual values:
		feasibility_iterate.update(penalty_iterate, step_size, combination_step);
		feasibility_iterate.update_dual(feasibility_step);
	}
	
	if (iter == maximum_iterations) {
		std::cout << "Failure - Did not converge" << std::endl;
	}
	std::cout << penalty_iterate.primal_values_ << std::endl;
	return 0;
}


